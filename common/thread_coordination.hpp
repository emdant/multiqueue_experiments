#ifndef THREAD_COORDINATION_HPP_7F7173D6
#define THREAD_COORDINATION_HPP_7F7173D6

#include "barrier.hpp"
#include "threading.hpp"

#include <pthread.h>

#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <ostream>
#include <utility>
#include <vector>

#ifndef L1_CACHE_LINESIZE
#define L1_CACHE_LINESIZE 64
#endif

namespace thread_coordination {

using duration_type = std::chrono::steady_clock::duration;

namespace detail {
template <class CharT, class Traits>
class ScopedOutput {
    using out_type = std::basic_ostream<CharT, Traits>;

   private:
    std::unique_lock<std::mutex> l_;
    out_type& out_;

   public:
    explicit ScopedOutput(std::mutex& mutex, out_type& out, int id) : out_{out} {
        l_ = std::unique_lock{mutex};
        out_ << "[Thread " << id << "] ";
    }

    template <typename T>
    friend out_type& operator<<(ScopedOutput const& so, T const& t) {
        so.out_ << t;
        return so.out_;
    }
};

struct SharedData {
    using clock_type = std::chrono::steady_clock;

    int num_threads;
    utils::barrier barrier;
    alignas(L1_CACHE_LINESIZE) std::atomic_size_t index{0};
    clock_type::time_point timestamp;
    std::mutex write_mutex;

    explicit SharedData(int n) : num_threads{n}, barrier{n} {
    }
};

}  // namespace detail

class Context {
    template <typename Affinity, typename Task>
    friend class TaskHandle;

    using duration_type = detail::SharedData::clock_type::duration;

    detail::SharedData& shared_data_;
    int id_;

    Context(detail::SharedData& sd, int id) : shared_data_{sd}, id_{id} {
    }

   public:
    [[nodiscard]] int get_num_threads() const noexcept {
        return shared_data_.num_threads;
    };

    [[nodiscard]] int get_id() const noexcept {
        return id_;
    }

    template <class CharT, class Traits>
    detail::ScopedOutput<CharT, Traits> write(std::basic_ostream<CharT, Traits>& out) const {
        return detail::ScopedOutput<CharT, Traits>(shared_data_.write_mutex, out, id_);
    }

    void synchronize() const {
        shared_data_.barrier.wait();
    }

    template <typename CompletionFunc>
    void synchronize(CompletionFunc&& f) const {
        shared_data_.barrier.wait(std::forward<CompletionFunc>(f));
    }

    template <typename Work, typename... Args>
    void execute_synchronized(Work work, Args&&... args) const {
        shared_data_.barrier.wait();
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait();
    }

    template <typename Work, typename... Args>
    void execute_synchronized_timed(duration_type& duration, Work work, Args&&... args) const {
        shared_data_.barrier.wait([this] { shared_data_.timestamp = std::chrono::steady_clock::now(); });
        work(std::forward<Args>(args)...);
        shared_data_.barrier.wait(
            [this, &duration] { duration = std::chrono::steady_clock::now() - shared_data_.timestamp; });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise(Iter begin, Iter end, Work work) const {
        using difference_type = typename std::iterator_traits<Iter>::difference_type;
        static constexpr difference_type block_size = 1 << 12;

        difference_type n = end - begin;
        shared_data_.barrier.wait();
        while (true) {
            auto block_begin =
                static_cast<difference_type>(shared_data_.index.fetch_add(block_size, std::memory_order_relaxed));
            if (block_begin >= n) {
                break;
            }
            auto block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait([this] { shared_data_.index.store(0, std::memory_order_relaxed); });
    }

    template <typename Iter, typename Work>
    void execute_synchronized_blockwise_timed(duration_type& duration, Iter begin, Iter end, Work work) const {
        using difference_type = typename std::iterator_traits<Iter>::difference_type;
        static constexpr difference_type block_size = 1 << 12;

        difference_type n = end - begin;
        shared_data_.barrier.wait([this] { shared_data_.timestamp = std::chrono::steady_clock::now(); });
        while (true) {
            auto block_begin =
                static_cast<difference_type>(shared_data_.index.fetch_add(block_size, std::memory_order_relaxed));
            if (block_begin >= n) {
                break;
            }
            auto block_end = std::min(block_begin + block_size, n);
            work(begin + block_begin, begin + block_end);
        }
        shared_data_.barrier.wait([this, &duration] {
            duration = std::chrono::steady_clock::now() - shared_data_.timestamp;
            shared_data_.index.store(0, std::memory_order_relaxed);
        });
    }
};

namespace affinity {

struct individual_cores {
    std::size_t stride = 1;
    std::size_t offset = 0;

    threading::thread_config operator()(int id) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(offset + static_cast<std::size_t>(id) * stride);
        return cfg;
    }
};

struct same_core {
    std::size_t core = 1;

    threading::thread_config operator()(int /*unused*/) const {
        threading::thread_config cfg;
        cfg.cpu_set.set(core);
        return cfg;
    }
};

}  // namespace affinity

template <typename Task, typename Affinity = affinity::individual_cores>
class TaskHandle {
    std::vector<threading::pthread> threads_;
    detail::SharedData shared_data_;

   public:
    template <typename... Args>
    explicit TaskHandle(Affinity affinity, int num_threads, Args... args)
        : threads_(static_cast<std::size_t>(num_threads)), shared_data_(num_threads) {
        for (std::size_t i = 0; i < threads_.size(); ++i) {
            Context ctx{shared_data_, static_cast<int>(i)};
            threads_[i] = threading::pthread(affinity(static_cast<int>(i)), Task::run, ctx, args...);
        }
    }

    template <typename... Args>
    explicit TaskHandle(int num_threads, Args... args) : TaskHandle(Affinity{}, num_threads, args...) {
    }

    TaskHandle(TaskHandle const&) = delete;
    TaskHandle(TaskHandle&&) noexcept = default;
    TaskHandle& operator=(TaskHandle const&) = delete;
    TaskHandle& operator=(TaskHandle&&) noexcept = default;
    ~TaskHandle() = default;

    void join() {
        for (auto& t : threads_) {
            t.join();
        }
    }
};

}  // namespace thread_coordination

#endif
