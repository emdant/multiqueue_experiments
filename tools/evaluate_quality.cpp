#include "cxxopts.hpp"
#include "replay_tree.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

struct InsertionLogEntry {
    std::uint64_t tick;
    unsigned long key;
    bool deleted;
};

struct DeletionLogEntry {
    unsigned int thread_id;
    std::uint64_t tick;
    bool failed;
    unsigned int insert_thread_id;
    unsigned long value;
};

std::istream& operator>>(std::istream& in, InsertionLogEntry& entry) {
    in >> entry.tick >> entry.key;
    entry.deleted = false;
    return in;
}

std::istream& operator>>(std::istream& in, DeletionLogEntry& entry) {
    in >> entry.thread_id >> entry.tick >> entry.insert_thread_id >>
        entry.value;
    entry.failed = false;
    return in;
}

std::ostream& operator<<(std::ostream& out, InsertionLogEntry const& entry) {
    out << entry.tick << ' ' << entry.key;
    return out;
}

std::ostream& operator<<(std::ostream& out, DeletionLogEntry const& entry) {
    out << entry.thread_id << ' ' << entry.tick << ' ' << entry.insert_thread_id
        << ' ' << entry.value;
    return out;
}

struct heap_entry {
    unsigned long key;
    unsigned int ins_thread_id;
    unsigned long elem_id;
};

struct extract_key {
    static unsigned long const& get(heap_entry const& e) { return e.key; }
};

int main(int argc, char* argv[]) {
    std::filesystem::path out_rank;
    std::filesystem::path out_delay;

    cxxopts::Options options(argv[0],
                             "Parses the logs generated by the generic test");
    options.positional_help("INSERT DELETION");
    // clang-format off
    options.add_options()
      ("v,verify", "Only verify the log")
      ("r,out-rank", "The output of the rank histogram", cxxopts::value<std::filesystem::path>(out_rank)->default_value("rank_histogram.txt"), "PATH")
      ("d,out-delay", "The output of the delay histogram", cxxopts::value<std::filesystem::path>(out_delay)->default_value("delay_histogram.txt"), "PATH")
      ("h,help", "Print this help");
    // clang-format on

    auto verify_only = false;
    try {
        auto result = options.parse(argc, argv);
        if (result.count("help")) {
            std::cerr << options.help() << std::endl;
            std::exit(0);
        }
        if (result.count("verify") > 0) {
            verify_only = true;
        }
    } catch (cxxopts::OptionParseException const& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }

    std::vector<std::vector<InsertionLogEntry>> insertions;
    std::vector<DeletionLogEntry> deletions;
    deletions.reserve(10'000'000);

    std::clog << "Reading quality log file from stdin..." << std::endl;
    unsigned int num_threads;
    std::cin >> num_threads;
    if (!std::cin || num_threads == 0) {
        std::cerr << "Line 1: "
                  << "Invalid number of threads" << std::endl;
        return 1;
    }
    insertions.resize(num_threads);

    char op;
    bool deleting = false;
    int line = 2;
    unsigned int num_corrections = 0;
    while (std::cin >> op) {
        if (op == 'i') {
            if (deleting) {
                std::cerr << "Line " << line << ": "
                          << "Insertion following a deletion" << std::endl;
                return 1;
            }
            unsigned int thread_id;
            std::cin >> thread_id;
            if (thread_id >= num_threads) {
                std::cerr << "Line " << line << ": "
                          << "Thread id " << thread_id
                          << " too high (Max: " << num_threads - 1 << ')'
                          << std::endl;
                return 1;
            }
            InsertionLogEntry entry;
            std::cin >> entry;
            if (!insertions[thread_id].empty() &&
                entry.tick < insertions[thread_id].back().tick) {
                std::cerr
                    << "Line " << line << ": "
                    << "Insertion\n\t" << entry
                    << "\nhappens before previous insertion of same thread"
                    << std::endl;
                return 1;
            }
            insertions[thread_id].push_back(entry);
        } else if (op == 'd') {
            deleting = true;
            DeletionLogEntry entry;
            std::cin >> entry;
            if (entry.thread_id >= num_threads) {
                std::cerr << "Line " << line << ": "
                          << "Thread id " << entry.thread_id
                          << " too high (Max: " << num_threads - 1 << ')'
                          << std::endl;
                return 1;
            }
            if (entry.insert_thread_id >= num_threads) {
                std::cerr << "Line " << line << ": "
                          << "Insert thread id " << entry.insert_thread_id
                          << " too high (Max: " << num_threads - 1 << ')'
                          << std::endl;
                return 1;
            }
            if (entry.value >= insertions[entry.insert_thread_id].size()) {
                std::cerr << "Line " << line << ": "
                          << "No insertion corresponding to deletion\n\t"
                          << entry << std::endl;
                return 1;
            }
            if (entry.tick <=
                insertions[entry.insert_thread_id][entry.value].tick) {
                ++num_corrections;
                /* insertions[entry.insert_thread_id][entry.value].tick =
                 * entry.tick - 1; */
                entry.tick =
                    insertions[entry.insert_thread_id][entry.value].tick + 1;
            }

            if (insertions[entry.insert_thread_id][entry.value].deleted) {
                std::cerr << "Line " << line << ": "
                          << "Insertion\n\t"
                          << insertions[entry.insert_thread_id][entry.value]
                          << "\n extracted twice" << std::endl;
                return 1;
            }
            deletions.push_back(entry);
            insertions[entry.insert_thread_id][entry.value].deleted = true;
        } else if (op == 'f') {
            deleting = true;
            DeletionLogEntry entry;
            std::cin >> entry.thread_id >> entry.tick;
            if (entry.thread_id >= num_threads) {
                std::cerr << "Line " << line << ": "
                          << "Thread id " << entry.thread_id
                          << " too high (Max: " << num_threads - 1 << ')'
                          << std::endl;
                return 1;
            }
            entry.failed = true;
            deletions.push_back(entry);
        } else {
            std::cerr << "Line " << line << ": "
                      << "Invalid operation: " << op << std::endl;
            return 1;
        }
        ++line;
    }
    std::clog << "Corrected insertion entries: " << num_corrections << '\n';
    if (verify_only) {
        std::clog << "Log is consistent" << std::endl;
        return 0;
    }

    std::clog << "Sorting operations..." << std::flush;
    /* std::for_each(insertions.begin(), insertions.end(), [](auto& list) { */
    /*   std::stable_sort(list.begin(), list.end(), [](auto const& lhs, auto
     * const& rhs) { */
    /*     return lhs.tick < rhs.tick; */
    /*   }); */
    /* }); */
    std::sort(
        deletions.begin(), deletions.end(),
        [](auto const& lhs, auto const& rhs) { return lhs.tick < rhs.tick; });
    std::clog << "done\n";
    std::vector<size_t> rank_histogram;
    std::vector<size_t> delay_histogram;

    std::clog << "Replaying operations...\n";
    ReplayTree<unsigned long, heap_entry, extract_key> replay_tree{};
    std::vector<size_t> insert_index(insertions.size(), 0);
    std::size_t failed_deletions = 0;

    for (auto it = deletions.begin(); it != deletions.end(); ++it) {
        // Inserting everything before next deletion
        for (unsigned int t = 0; t < insertions.size(); ++t) {
            while (insert_index[t] < insertions[t].size() &&
                   insertions[t][insert_index[t]].tick < it->tick) {
                replay_tree.insert(
                    {insertions[t][insert_index[t]].key, t,
                     static_cast<unsigned long>(insert_index[t])});
                ++insert_index[t];
            }
        }

        if (it->failed) {
            if (!replay_tree.empty()) {
                ++failed_deletions;
                size_t rank_error = replay_tree.size();
                if (rank_histogram.size() <= rank_error) {
                    rank_histogram.resize(rank_error + 1, 0);
                }
                ++rank_histogram[rank_error];
                replay_tree.increase_global_delay();
            }
            continue;
        }

        auto key = insertions[it->insert_thread_id][it->value].key;
        auto [start, end] = replay_tree.equal_range(key);
        while (start != end && (start->ins_thread_id != it->insert_thread_id ||
                                start->elem_id != it->value)) {
            ++start;
        }
        if (start == end) {
            std::cerr << "Element\n\t" << *it
                      << "\nis not in the heap at deletion time" << std::endl;
            return 1;
        }
        size_t rank_error = replay_tree.get_rank(key);
        if (rank_histogram.size() <= rank_error) {
            rank_histogram.resize(rank_error + 1, 0);
        }
        ++rank_histogram[rank_error];
        auto [success, delay] = replay_tree.erase(start);
        if (!success || delay < 0) {
            std::cerr << delay << '\n';
            std::cerr << "Inconsistent replay tree state" << std::endl;
            return 1;
        }
        if (delay_histogram.size() <= static_cast<std::size_t>(delay)) {
            delay_histogram.resize(static_cast<std::size_t>(delay) + 1, 0);
        }
        ++delay_histogram[static_cast<std::size_t>(delay)];
        replay_tree.increase_delay(key);
        if (static_cast<std::size_t>(std::distance(deletions.begin(), it)) %
                (deletions.size() / 100) ==
            0) {
            std::clog << "\rProcessed " << std::setprecision(3)
                      << 100.0 *
                             static_cast<double>(
                                 std::distance(deletions.begin(), it)) /
                             static_cast<double>(deletions.size())
                      << "%";
        }
    }
    std::clog << "\rProcessing done         " << std::endl;
    std::clog << "Failed deletions: " << failed_deletions << std::endl;
    std::clog << "Writing histograms..." << std::flush;
    {
        auto out_f = std::ofstream{out_rank};
        for (size_t i = 0; i < rank_histogram.size(); ++i) {
            if (rank_histogram[i] > 0) {
                out_f << i << " " << rank_histogram[i] << '\n';
            }
        }
        out_f.close();
    }
    {
        auto out_f = std::ofstream{out_delay};
        for (size_t i = 0; i < delay_histogram.size(); ++i) {
            if (delay_histogram[i] > 0) {
                out_f << i << " " << delay_histogram[i] << '\n';
            }
        }
        out_f.close();
    }
    std::clog << "done" << std::endl;
    return 0;
}
