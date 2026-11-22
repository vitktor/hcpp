#pragma once
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace hc
{

// Number of worker threads to use for a requested count: 0 means one per
// hardware thread (at least 1).
inline int resolve_threads(int requested)
{
  if (requested > 0)
    return requested;
  return std::max(1u, std::thread::hardware_concurrency());
}

// Calls f(thread_id, i) for every i in [0, n), spread over at most
// num_threads threads (0 = one per hardware thread). thread_id is in
// [0, min(num_threads, n)) and identifies the calling worker, so f can use
// per-thread state (e.g. a homotopy that is mutated while tracking).
//
// Indices are handed out one at a time from a shared atomic counter rather
// than in fixed blocks: path tracking times vary by orders of magnitude, so
// static partitioning would leave threads idle behind a few hard paths.
//
// f must be safe to call concurrently for distinct i; writing to a
// preallocated output slot per index needs no locking. If any call throws,
// the remaining indices are skipped and the first exception is rethrown
// after all workers have joined. With one thread (or n <= 1), f runs inline
// on the calling thread.
template <typename F>
void parallel_for(std::size_t n, int num_threads, F&& f)
{
  std::size_t workers = std::min<std::size_t>(resolve_threads(num_threads), n);
  if (workers <= 1)
  {
    for (std::size_t i = 0; i < n; ++i)
      f(0, i);
    return;
  }

  std::atomic<std::size_t> next{0};
  std::atomic<bool> failed{false};
  std::vector<std::exception_ptr> errors(workers);

  auto work = [&](int tid) {
    try
    {
      for (std::size_t i = next.fetch_add(1, std::memory_order_relaxed);
           i < n && !failed.load(std::memory_order_relaxed);
           i = next.fetch_add(1, std::memory_order_relaxed))
        f(tid, i);
    }
    catch (...)
    {
      errors[tid] = std::current_exception();
      failed.store(true, std::memory_order_relaxed);
    }
  };

  // std::thread rather than std::jthread: Apple's libc++ only ships jthread
  // behind -fexperimental-library.
  std::vector<std::thread> threads;
  threads.reserve(workers - 1);
  try
  {
    for (std::size_t t = 1; t < workers; ++t)
      threads.emplace_back(work, static_cast<int>(t));
  }
  catch (...)
  {
    // Thread creation failed: stop the ones already running before
    // propagating, as destroying a joinable std::thread terminates.
    failed.store(true);
    for (auto& th : threads)
      th.join();
    throw;
  }
  work(0); // the calling thread is worker 0
  for (auto& th : threads)
    th.join();

  for (auto& e : errors)
    if (e)
      std::rethrow_exception(e);
}

} // namespace hc
