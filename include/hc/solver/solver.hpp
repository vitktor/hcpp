#pragma once

#include <hc/core/system.hpp>
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
template <typename Scalar>
class Solver
{
public:
  explicit Solver(SolverOptions<Scalar> options = {}) : options_(options) {}

  Result<Scalar> solve(const System<Scalar>& target) const
  {
    TotalDegreeSystem<Scalar> start(target);
    StraightLineHomotopy<Scalar> H(start.system(), target);
    EulerPredictor<Scalar> predictor;
    Tracker<Scalar> tracker(H, predictor, options_.corrector_tol, options_.corrector_max_iters,
                            options_.dt_init, options_.dt_min, options_.dt_max,
                            options_.max_steps, options_.max_norm);

    Result<Scalar> result;
    result.paths.reserve(start.start_solutions().size());
    for (const auto& x0 : start.start_solutions())
      result.paths.push_back(tracker.track(x0, Scalar(0), Scalar(1)));
    return result;
  }

private:
  SolverOptions<Scalar> options_;
};

} // namespace hc
