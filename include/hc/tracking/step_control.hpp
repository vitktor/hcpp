#pragma once
#include <algorithm>

namespace hc
{

// Tracks a step size magnitude: grows it after enough consecutive predict+correct
// successes, shrinks it on any failure, always clamped to [dt_min, dt_max].
// Direction (increasing/decreasing t) is the caller's concern, not this class's.
template <typename Scalar>
class AdaptiveStepController
{
public:
  AdaptiveStepController(Scalar dt_init, Scalar dt_min, Scalar dt_max,
                         int success_streak_for_growth = 3,
                         Scalar growth_factor = Scalar(2), Scalar shrink_factor = Scalar(0.5))
      : dt_(dt_init), dt_min_(dt_min), dt_max_(dt_max),
        success_streak_for_growth_(success_streak_for_growth),
        growth_factor_(growth_factor), shrink_factor_(shrink_factor)
  {
  }

  Scalar step_size() const { return dt_; }

  Scalar on_success()
  {
    if (++success_streak_ >= success_streak_for_growth_)
    {
      dt_ = std::min(dt_ * growth_factor_, dt_max_);
      success_streak_ = 0;
    }
    return dt_;
  }

  Scalar on_failure()
  {
    success_streak_ = 0;
    dt_ = std::max(dt_ * shrink_factor_, dt_min_);
    return dt_;
  }

private:
  Scalar dt_;
  Scalar dt_min_;
  Scalar dt_max_;
  int success_streak_for_growth_;
  Scalar growth_factor_;
  Scalar shrink_factor_;
  int success_streak_ = 0;
};

} // namespace hc
