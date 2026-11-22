#pragma once
#include <cstddef>
#include <mutex>

#include <hc/core/system.hpp>
#include <hc/core/types.hpp>
#include <hc/parallel.hpp>
#include <hc/homotopy/straight_line.hpp>
#include <hc/solver/result.hpp>
#include <hc/solver/solver_options.hpp>
#include <hc/start_systems/total_degree.hpp>
#include <hc/tracking/runge_kutta.hpp>
#include <hc/tracking/tracker.hpp>

namespace hc
{

// Solves target=0 by building a TotalDegreeSystem start system, tracking
// every one of its start solutions to target via a StraightLineHomotopy
// (its default random gamma providing the genericity needed to avoid path
// crossings/singularities), and collecting the results.
//
// Paths are tracked in parallel on options.threads threads. The homotopy,
// predictor and tracker are shared read-only (all their tracking methods are
// const), and each path writes only its own slot of result.paths, so
// results come back in start-solution order regardless of thread count.
//
// With options.verbose, each worker reports its finished path to
// options.log. Those writes share one stream and one counter, so they go
// through a mutex: unguarded, concurrent writes to a std::ostream are a data
// race and lines from different threads would interleave.
template <typename Scalar>
class Solver
{
public:
  explicit Solver(SolverOptions<Scalar> options = {}) : options_(options) {}

  Result<Scalar> solve(const System<Scalar>& target) const
  {
    TotalDegreeSystem<Scalar> start(target);
    StraightLineHomotopyOptions<Scalar> h_options;
    if (options_.gamma)
      h_options.gamma = *options_.gamma;
    StraightLineHomotopy<Scalar> H(start.system(), target, h_options);
    EulerPredictor<Scalar> predictor;
    Tracker<Scalar> tracker(H, predictor, options_.corrector_tol, options_.corrector_max_iters,
                            options_.dt_init, options_.dt_min, options_.dt_max,
                            options_.max_steps, options_.max_norm);

    const auto& starts = start.start_solutions();
    Result<Scalar> result;
    result.paths.resize(starts.size());
    std::mutex log_mutex;
    std::size_t finished = 0; // guarded by log_mutex
    parallel_for(starts.size(), options_.threads, [&](int, size_t i) {
      result.paths[i] = tracker.track(starts[i], Scalar(0), Scalar(1));
      if (options_.verbose)
      {
        const auto& path = result.paths[i];
        // Increment and print under one lock, so the counter in the output
        // runs 1, 2, ..., n in order.
        std::lock_guard lock(log_mutex);
        ++finished;
        *options_.log << "[" << finished << "/" << starts.size() << "] path " << i << ": "
                      << to_string(path.status) << ", " << path.steps << " steps\n";
      }
    });
    return result;
  }

private:
  SolverOptions<Scalar> options_;
};

} // namespace hc
