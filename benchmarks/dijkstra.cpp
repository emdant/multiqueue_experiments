#include "util/benchmark.h"
#include "util/build_info.hpp"
#include "util/selector.hpp"
#include "util/termination_detection.hpp"
#include "util/thread_coordination.hpp"

#include <cxxopts.hpp>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <x86intrin.h>
#include <array>
#include <atomic>
#include <cassert>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <random>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using pq_type = PQ<true, unsigned long, unsigned long>;
using handle_type = pq_type::handle_type;
using node_type = pq_type::value_type;

struct Settings {
    int num_threads = 4;
    int trials = 1;
    int delta = 1;
    std::filesystem::path graph_file;
    unsigned int seed = 1;
    pq_type::settings_type pq_settings{};
};

void register_cmd_options(Settings& settings, cxxopts::Options& cmd) {
    // clang-format off
    cmd.add_options()
        ("j,threads", "The number of threads", cxxopts::value<int>(settings.num_threads), "NUMBER")
        ("n,trials", "The number of trials", cxxopts::value<int>(settings.trials), "NUMBER")
        ("d,delta", "Delta", cxxopts::value<int>(settings.delta), "NUMBER")
        ("graph", "The input graph", cxxopts::value<std::filesystem::path>(settings.graph_file), "PATH");
    // clang-format on
    settings.pq_settings.register_cmd_options(cmd);
    cmd.parse_positional({"graph"});
}

void write_settings_human_readable(Settings const& settings, std::ostream& out) {
    out << "Threads: " << settings.num_threads << '\n';
    out << "Trials: " << settings.trials << '\n';
    out << "Delta: " << settings.delta << '\n';
    out << "Graph: " << settings.graph_file << '\n';
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
    CumulativeTimer push_timer;
    CumulativeTimer pop_timer;
    // long long pushed_nodes{0};
    // long long ignored_nodes{0};
    // long long processed_nodes{0};
    // long long visits{0};
};

struct alignas(L1_CACHE_LINE_SIZE) AtomicDistance {
    std::atomic<long long> value{std::numeric_limits<long long>::max()};
};

struct SharedData {
    WGraph& graph;
    std::vector<AtomicDistance> distances;
    termination_detection::TerminationDetection termination_detection;
};

// void process_node(node_type const& node, handle_type& handle, Counter& counter, SharedData& data) {
//     auto current_distance = data.distances[node.second].value.load(std::memory_order_relaxed);
//     if (static_cast<long long>(node.first) > current_distance) {
//         // ++counter.ignored_nodes;
//         return;
//     }
//     for (WNode wn : data.graph.out_neigh(node.second)) {
//         // ++counter.visits;
//         auto target = wn.v;
//         auto d = current_distance + wn.w;
//         auto old_d = data.distances[target].value.load(std::memory_order_relaxed);
//         while (d < old_d) {
//             if (data.distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_relaxed)) {
//                 counter.push_timer.Start();
//                 handle.push({d, target});
//                 counter.push_timer.Stop();
//                 // ++counter.pushed_nodes;
//                 break;
//             }
//         }
//     }
//     // ++counter.processed_nodes;
// }

void process_node(node_type const& node, handle_type& handle, Counter& counter, SharedData& data, int delta) {
    auto current_distance = data.distances[node.second].value.load(std::memory_order_relaxed);
    if (static_cast<long long>(node.first) * delta > current_distance) {
        // ++counter.ignored_nodes;
        return;
    }
    for (WNode wn : data.graph.out_neigh(node.second)) {
        auto target = wn.v;
        auto d = current_distance + wn.w;
        auto old_d = data.distances[target].value.load(std::memory_order_relaxed);
        while (d < old_d) {
            if (data.distances[target].value.compare_exchange_weak(old_d, d, std::memory_order_relaxed)) {
                // counter.push_timer.Start();
                handle.push({d / delta, target});
                // counter.push_timer.Stop();
                // ++counter.pushed_nodes;
                break;
            }
        }
    }
    // ++counter.processed_nodes;
}

[[gnu::noinline]] Counter benchmark_thread(thread_coordination::Context& thread_context, pq_type& pq, SharedData& data,
                                           NodeID src, int delta) {
    Counter counter;
    auto handle = pq.get_handle();
    if (thread_context.id() == 0) {
        data.distances[src].value = 0;
        handle.push({0, src});
        // ++counter.pushed_nodes;
    }
    thread_context.synchronize();
    std::optional<node_type> node;
    while (data.termination_detection.repeat([&]() {
        // counter.pop_timer.Start();
        node = handle.try_pop();
        // counter.pop_timer.Stop();
        return node.has_value();
    })) {
        process_node(*node, handle, counter, data, delta);
    }
    thread_context.synchronize();
    return counter;
}

void run_benchmark(WGraph& g, Settings const& settings) {
    SourcePicker<WGraph> sp(g);
    NodeID src = sp.PickNext();

    SharedData shared_data{g, {}, termination_detection::TerminationDetection(settings.num_threads)};
    std::cout << "\nSource: " << src << std::endl;

    shared_data.distances = std::vector<AtomicDistance>(shared_data.graph.num_nodes());

    std::vector<Counter> thread_counter(static_cast<std::size_t>(settings.num_threads));
    auto pq = pq_type(settings.num_threads, shared_data.graph.num_nodes(), settings.pq_settings);
    auto start_time = std::chrono::steady_clock::now();
    thread_coordination::Dispatcher dispatcher{settings.num_threads, [&](auto ctx) {
                                                   auto t_id = static_cast<std::size_t>(ctx.id());
                                                   thread_counter[t_id] =
                                                       benchmark_thread(ctx, pq, shared_data, src, settings.delta);
                                               }};
    dispatcher.wait();
    auto end_time = std::chrono::steady_clock::now();

    // double push_time = 0, pop_time = 0;
    // for (auto i = 0; i < settings.num_threads; i++) {
    //     push_time += thread_counter[i].push_timer.Seconds();
    //     pop_time += thread_counter[i].pop_timer.Seconds();
    // }
    // push_time /= settings.num_threads;
    // pop_time /= settings.num_threads;

    // auto total_counts =
    //     std::accumulate(thread_counter.begin(), thread_counter.end(), Counter{}, [](auto sum, auto const&
    //     counter) {
    //         sum.pushed_nodes += counter.pushed_nodes;
    //         sum.processed_nodes += counter.processed_nodes;
    //         sum.ignored_nodes += counter.ignored_nodes;
    //         sum.visits += counter.visits;
    //         return sum;
    //     });
    auto longest_distance =
        std::max_element(shared_data.distances.begin(), shared_data.distances.end(), [](auto const& a, auto const& b) {
            auto a_val = a.value.load(std::memory_order_relaxed);
            auto b_val = b.value.load(std::memory_order_relaxed);
            if (b_val == std::numeric_limits<long long>::max()) {
                return false;
            }
            if (a_val == std::numeric_limits<long long>::max()) {
                return true;
            }
            return a_val < b_val;
        })->value.load();
    std::clog << "= Results =\n";
    std::clog << "Time (s): " << std::fixed << std::setprecision(6)
              << std::chrono::duration<double>(end_time - start_time).count() << '\n';
    std::clog << "Longest distance: " << longest_distance << '\n';
    // std::clog << "Push time: " << push_time << '\n';
    // std::clog << "Pop time: " << pop_time << '\n';
    // std::clog << "Processed nodes: " << total_counts.processed_nodes << '\n';
    // std::clog << "Ignored nodes: " << total_counts.ignored_nodes << '\n';
    // std::clog << "Number of relaxations: " << total_counts.visits << '\n';

    // if (total_counts.processed_nodes + total_counts.ignored_nodes != total_counts.pushed_nodes) {
    //     std::cerr << "Warning: Not all nodes were popped\n";
    //     std::cerr << "Probably the priority queue discards duplicate keys\n";
    // }
    // std::cout << '{';
    // std::cout << std::quoted("settings") << ':';
    // write_settings_json(settings, std::cout);
    // std::cout << ',';
    // std::cout << std::quoted("graph") << ':';
    // std::cout << '{';
    // std::cout << std::quoted("num_nodes") << ':' << shared_data.graph.num_nodes() << ',';
    // std::cout << std::quoted("num_edges") << ':' << shared_data.graph.num_edges();
    // std::cout << '}' << ',';
    // std::cout << std::quoted("results") << ':';
    // std::cout << '{';
    // std::cout << std::quoted("time_ns") << ':' << std::chrono::nanoseconds{end_time - start_time}.count() << ',';
    // std::cout << std::quoted("longest_distance") << ':' << longest_distance << ',';
    // std::cout << std::quoted("processed_nodes") << ':' << total_counts.processed_nodes << ',';
    // std::cout << std::quoted("ignored_nodes") << ':' << total_counts.ignored_nodes;
    // std::cout << '}';
    // std::cout << '}' << '\n';
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

    std::clog << "= Running benchmark =\n";
    for (auto i = 0; i < settings.trials; i++)
        run_benchmark(g, settings);
    return EXIT_SUCCESS;
}
