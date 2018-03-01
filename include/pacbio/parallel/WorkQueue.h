// Author: Lance Hepler

#pragma once

#include <condition_variable>
#include <exception>
#include <future>
#include <mutex>
#include <queue>

#include <boost/optional.hpp>

namespace PacBio {
namespace Parallel {

template <typename T>
class WorkQueue
{
private:
    typedef boost::optional<std::packaged_task<T(void)>> TTask;
    typedef boost::optional<std::future<T>> TFuture;

public:
    WorkQueue(const size_t size) : exc{nullptr}, sz{size}
    {
        for (size_t i = 0; i < size; ++i) {
            threads.emplace_back(std::thread([this]() {
                try {
                    while (auto task = PopTask()) {
                        (*task)();
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> g(m);
                        exc = std::current_exception();
                    }
                    popped.notify_all();
                }
            }));
        }
    }

    ~WorkQueue()
    {
        for (auto& thread : threads) {
            thread.join();
        }

        // wait to see if there's a final exception, throw if so..
        {
            std::lock_guard<std::mutex> g(m);
            if (exc) std::rethrow_exception(exc);
        }
    }

    template <typename F, typename... Args>
    void ProduceWith(F&& f, Args&&... args)
    {
        std::packaged_task<T(void)> task{
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)};

        {
            std::unique_lock<std::mutex> lk(m);
            popped.wait(lk, [&task, this]() {
                if (exc) std::rethrow_exception(exc);

                if (head.size() < sz) {
                    head.emplace(std::move(task));
                    return true;
                }

                return false;
            });
        }
        pushed.notify_all();
    }

    template <typename F, typename... Args>
    bool ConsumeWith(F&& cont, Args&&... args)
    {
        TFuture fut(boost::none);

        {
            std::unique_lock<std::mutex> lk(m);
            popped.wait(lk, [&fut, this]() {
                if (tail.empty()) return false;

                if ((fut = std::move(tail.front()))) {
                    tail.pop();
                }

                return true;
            });
        }

        try {
            if (!fut) return false;

            cont(std::forward<Args>(args)..., std::move(fut->get()));
            return true;
        } catch (...) {
            {
                std::lock_guard<std::mutex> g(m);
                exc = std::current_exception();
            }
            popped.notify_all();
        }
        return false;
    }

    void Finalize()
    {
        {
            std::lock_guard<std::mutex> g(m);
            head.emplace(boost::none);
        }
        pushed.notify_all();
    }

private:
    TTask PopTask()
    {
        TTask task(boost::none);

        {
            std::unique_lock<std::mutex> lk(m);
            pushed.wait(lk, [&task, this]() {
                if (head.empty()) return false;

                if ((task = std::move(head.front()))) {
                    head.pop();
                    tail.emplace(task->get_future());
                } else
                    tail.emplace(boost::none);

                return true;
            });
        }
        popped.notify_all();

        return task;
    }

    std::vector<std::thread> threads;
    std::queue<TTask> head;
    std::queue<TFuture> tail;
    std::condition_variable popped;
    std::condition_variable pushed;
    std::exception_ptr exc;
    std::mutex m;
    size_t sz;
};

}  // namespace parallel
}  // namespace PacBio
