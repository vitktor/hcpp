#pragma once
#include <algorithm>
#include <complex>
#include <limits>
#include <vector>

#include <hc/core/types.hpp>
#include <hc/homotopy/homotopy.hpp>
#include <hc/tracking/corrector.hpp>
#include <hc/tracking/predictor.hpp>
#include <hc/tracking/step_control.hpp>

namespace hc
{

template <typename Scalar>
struct TrackResult
{
  std::vector<std::complex<Scalar>> solution;
  bool success;
  PathStatus status;
  Scalar t;
  int steps;
};

// Tracks a single path of H(x,t)=0 from t_start to t_end via a predict/correct
// loop with adaptive step size. Stateless across calls: each track() call builds
// its own corrector and step controller, so one Tracker can track many
// independent paths (e.g. from multiple start solutions) without interference.
template <typename Scalar>
class Tracker
{
public:
  Tracker(const Homotopy<Scalar>& H, const Predictor<Scalar>& predictor,
         Scalar corrector_tol, int corrector_max_iters,
         Scalar dt_init, Scalar dt_min, Scalar dt_max, int max_steps,
         Scalar max_norm = Scalar(1e6))
      : H_(H), predictor_(predictor), corrector_tol_(corrector_tol),
        corrector_max_iters_(corrector_max_iters), dt_init_(dt_init),
        dt_min_(dt_min), dt_max_(dt_max), max_steps_(max_steps), max_norm_(max_norm)
  {
  }

  TrackResult<Scalar> track(std::vector<std::complex<Scalar>> x, Scalar t_start,
                           Scalar t_end) const
  {
    NewtonCorrector<Scalar> corrector(corrector_tol_, corrector_max_iters_);
    AdaptiveStepController<Scalar> step(dt_init_, dt_min_, dt_max_);
    Scalar direction = t_end > t_start ? Scalar(1) : Scalar(-1);
    Scalar t = t_start;
    Scalar eps = std::numeric_limits<Scalar>::epsilon() * Scalar(100);

    for (int i = 0; i < max_steps_; ++i)
    {
      Scalar remaining = direction * (t_end - t);
      if (remaining <= eps)
        return {std::move(x), true, PathStatus::Success, t, i};

      Scalar dt = direction * std::min(step.step_size(), remaining);

      std::vector<std::complex<Scalar>> predicted;
      predictor_.predict(x, t, dt, H_, predicted);
      Scalar t_new = t + dt;

      auto result = corrector.correct(predicted, H_, t_new);
      if (result.converged)
      {
        x = std::move(predicted);
        t = t_new;
        step.on_success();

        if (diverged(x))
          return {std::move(x), false, PathStatus::Diverged, t, i + 1};
      }
      else
      {
        step.on_failure();
        if (step.step_size() <= dt_min_)
          return {std::move(x), false, PathStatus::Failed, t, i + 1};
      }
    }

    return {std::move(x), false, PathStatus::Failed, t, max_steps_};
  }

private:
  bool diverged(const std::vector<std::complex<Scalar>>& x) const
  {
    for (const auto& xi : x)
      if (std::abs(xi) > max_norm_)
        return true;
    return false;
  }

  const Homotopy<Scalar>& H_;
  const Predictor<Scalar>& predictor_;
  Scalar corrector_tol_;
  int corrector_max_iters_;
  Scalar dt_init_;
  Scalar dt_min_;
  Scalar dt_max_;
  int max_steps_;
  Scalar max_norm_;
};

} // namespace hc
