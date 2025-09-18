#include "util/benchmark.h"
#include "util/build_info.hpp"
#include "util/selector.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"

#include <algorithm>
#include <cstdint>
#include <cxxopts.hpp>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <x86intrin.h>
#include <atomic>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#ifdef USE_FLOAT
typedef float WeightT;
#else
typedef int32_t WeightT;
#endif

using pq_type = PQ<true, WeightT, NodeID>;
using handle_type = pq_type::handle_type;
using node_type = pq_type::value_type;

struct Settings {
    NodeID src;
    int num_threads = 4;
    int sources = 1;
    int trials = 1;
    std::filesystem::path graph_file;
    std::filesystem::path sources_file = "";
    unsigned int seed = 1;
    pq_type::settings_type pq_settings{};
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    // clang-format off
    cmd.add_options()
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("S,sources", "The number of source nodes", cxxopts::value<int>(settings.sources), "NUMBER")
        ("n,trials", "The number of trials per source", cxxopts::value<int>(settings.trials), "NUMBER")
        ("graph", "The input graph", cxxopts::value<std::filesystem::path>(settings.graph_file), "PATH")
        ("z,sources_file", "The input sources", cxxopts::value<std::filesystem::path>(settings.sources_file), "PATH");
    // clang-format on
    settings.pq_settings.register_cmd_options(cmd);
    cmd.parse_positional({"graph"});
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Trials: " << settings.trials << '\n';
    out << "Graph: " << settings.graph_file << '\n';
    out << "Sources: " << ((settings.sources_file.string() == "") ? "randomly generated" : settings.sources_file)
        << '\n';
    settings.pq_settings.write_human_readable(out);
}

void write_settings_json(Settings const& settings, std::ostream& out) {
    out << '{';
    out << std::quoted("num_threads") << ':' << settings.num_threads << ',';
    out << std::quoted("num_trials") << ':' << settings.trials << ',';
    out << std::quoted("graph_file") << ':' << settings.graph_file << ',';
    out << std::quoted("seed") << ':' << settings.seed << ',';
    out << std::quoted("pq") << ':';
    settings.pq_settings.write_json(out);
    out << '}';
}

struct Counter {
#ifdef COUNT_TIME
    CumulativeTimer push_timer;
    CumulativeTimer pop_timer;
#endif
};

struct alignas(L1_CACHE_LINE_SIZE) AtomicDistance {
    std::atomic<WeightT> value{std::numeric_limits<WeightT>::max()};
};

struct SharedData {
    WGraph& graph;
    std::vector<AtomicDistance> distances;
    termination_detection::TerminationDetection termination_detection;
};

void process_node(node_type const& node, handle_type& handle, Counter& counter, SharedData& data) {
    WeightT current_distance = data.distances[node.second].value.load(std::memory_order_relaxed);
    if (static_cast<WeightT>(node.first) > current_distance) {
        return;
    }
    for (WNode wn : data.graph.out_neigh(node.second)) {
        auto target = wn.v;
        auto d = static_cast<WeightT>(node.first) + wn.w;
        auto old_d = data.distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (data.distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_relaxed)) {
#ifdef COUNT_TIME
                counter.push_timer.Start();
#endif
                handle.push({d, target});
#ifdef COUNT_TIME
                counter.push_timer.Stop();
#endif
                break;
            }
        }
    }
}

[[gnu::noinline]] Counter benchmark_thread(thread_coordination::Context& thread_context, pq_type& pq, SharedData& data,
                                           NodeID src) {
    Counter counter;
    auto handle = pq.get_handle();
    if (thread_context.id() == 0) {
        data.distances[src].value = 0;
        handle.push({0, src});
    }
    thread_context.synchronize();
    std::optional<node_type> node;
    while (data.termination_detection.repeat([&]() {
#ifdef COUNT_TIME
        counter.pop_timer.Start();
#endif
        node = handle.try_pop();
#ifdef COUNT_TIME
        counter.pop_timer.Stop();
#endif
        return node.has_value();
    })) {
        process_node(*node, handle, counter, data);
    }
    thread_context.synchronize();
    return counter;
}

void run_benchmark(WGraph& g, Settings const& settings) {
    SharedData shared_data{g, {}, termination_detection::TerminationDetection(settings.num_threads)};

    shared_data.distances = std::vector<AtomicDistance>(shared_data.graph.num_nodes());

    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.num_threads));
    auto pq = pq_type(settings.num_threads, shared_data.graph.num_nodes(), settings.pq_settings);
    auto start_time = std::chrono::steady_clock::now();
    thread_coordination::Dispatcher dispatcher{settings.num_threads, [&](auto ctx) {
                                                   auto t_id = static_cast<std::size_t>(ctx.id());
                                                   thread_counter[t_id] =
                                                       benchmark_thread(ctx, pq, shared_data, settings.src);
                                               }};
    dispatcher.wait();
    auto end_time = std::chrono::steady_clock::now();
#ifdef COUNT_TIME
    double push_time = 0, pop_time = 0;
    for (auto i = 0; i < settings.num_threads; i++) {
        push_time += thread_counter[i].push_timer.Seconds();
        pop_time += thread_counter[i].pop_timer.Seconds();
    }
    push_time /= settings.num_threads;
    pop_time /= settings.num_threads;
#endif
    auto longest_distance =
        std::max_element(shared_data.distances.begin(), shared_data.distances.end(), [](auto const& a, auto const& b) {
            auto a_val = a.value.load(std::memory_order_relaxed);
            auto b_val = b.value.load(std::memory_order_relaxed);
            if (b_val == std::numeric_limits<WeightT>::max()) {
                return false;
            }
            if (a_val == std::numeric_limits<WeightT>::max()) {
                return true;
            }
            return a_val < b_val;
        })->value.load();
    std::clog << "Time (s): " << std::fixed << std::setprecision(6)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    NodeID num_reached = std::count_if(shared_data.distances.begin(), shared_data.distances.end(), [&](auto const& d) {
        return d.value.load(std::memory_order_relaxed) != std::numeric_limits<WeightT>::max();
    });
    std::clog << "Nodes reached: " << num_reached << '\n';
    std::clog << "Longest distance: " << longest_distance << std::endl;

#ifdef COUNT_TIME
    std::clog << "Push time: " << push_time << '\n';
    std::clog << "Pop time: " << pop_time << '\n';
#endif
}

int main(int argc, char* argv[]) {
    write_build_info(std::clog);
    std::clog << '\n';

    std::clog << "= Priority queue =\n";
    pq_type::write_human_readable(std::clog);
    std::clog << '\n';

    std::clog << "= Command line =\n";
    for (int i = 0; i < argc; ++i) {
        std::clog << argv[i];
        if (i != argc - 1) {
            std::clog << ' ';
        }
    }
    std::clog << '\n' << '\n';

    cxxopts::Options cmd(argv[0]);
    cmd.add_options()("h,help", "Print this help");
    Settings settings{};
    register_cmd_options(settings, cmd);

    try {
        auto args = cmd.parse(argc, argv);
        if (args.count("help") > 0) {
            std::cerr << cmd.help() << '\n';
            return EXIT_SUCCESS;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << "Error parsing command line: " << e.what() << '\n';
        std::cerr << "Use --help for usage information" << '\n';
        return EXIT_FAILURE;
    }

    std::clog << "= Settings =\n";
    write_settings_human_readable(settings, std::clog);
    std::clog << '\n';

    std::clog << "Reading graph...\n";
    WGraph g = MakeWeightedGraph(settings.graph_file.string());
    g.PrintStats();

    SourcePicker<WGraph> sp(g, settings.sources_file.string());
    std::clog << "= Running benchmark =\n";
    for (auto s = 0; s < settings.sources; s++) {
        settings.src = sp.PickNext();
        std::cout << "\nSource: " << settings.src << std::endl;

        for (auto i = 0; i < settings.trials; i++) {
            run_benchmark(g, settings);
        }
    }
    return EXIT_SUCCESS;
}
