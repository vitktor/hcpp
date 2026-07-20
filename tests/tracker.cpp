#include <cmath>

#include <gtest/gtest.h>
#include <hc/homotopy/straight_line.hpp>
#include <hc/tracking/runge_kutta.hpp>
#include <hc/tracking/tracker.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(Tracker, TracksNonlinearPathFromStartToTarget) {
  // G = x^2 - 4 (root 2, the known start solution at t=0)
  // F = x^2 - 9 (root 3, the target solution at t=1)
  // H(x,t) = x^2 - (4 + 5t); track t: 0 -> 1
  Variable x("x");
  System<double> start({pow(x, 2) - 4.0}, {x});
  System<double> target({pow(x, 2) - 9.0}, {x});
  StraightLineHomotopy<double> H(start, target);
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor, /*corrector_tol=*/1e-12, /*corrector_max_iters=*/20,
                          /*dt_init=*/0.05, /*dt_min=*/1e-8, /*dt_max=*/0.1,
                          /*max_steps=*/200);

  auto result = tracker.track({cd(2.0, 0.0)}, 0.0, 1.0);

  EXPECT_TRUE(result.success);
  EXPECT_NEAR(result.t, 1.0, 1e-9);
  EXPECT_NEAR(result.solution[0].real(), 3.0, 1e-8);
  EXPECT_NEAR(result.solution[0].imag(), 0.0, 1e-8);
}

TEST(Tracker, TracksCircleAndLineFromStartToTarget) {
  // G = { x^2 - 1, y^2 - 1 } (known start solution (1,1) at t=0)
  // F = { x^2 + y^2 - 1, x - y } (target solution (1/sqrt2, 1/sqrt2) at t=1)
  // Without a gamma trick this path hits a non-generic singularity around
  // t=1/3 and fails (see examples/circle_and_line); gamma=0.6+0.8i steers
  // it off the real axis and it converges cleanly in 12 steps.
  Variable x("x"), y("y");
  auto g1 = pow(x, 2) - 1.0;
  auto g2 = pow(y, 2) - 1.0;
  auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
  auto f2 = x - y;
  System<double> start({g1, g2}, {x, y});
  System<double> target({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(start, target, cd(0.6, 0.8));
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor, 1e-12, 20, 0.05, 1e-8, 0.1, 200);

  auto result = tracker.track({cd(1.0, 0.0), cd(1.0, 0.0)}, 0.0, 1.0);

  double expected = 1.0 / std::sqrt(2.0);
  EXPECT_TRUE(result.success);
  EXPECT_NEAR(result.solution[0].real(), expected, 1e-8);
  EXPECT_NEAR(result.solution[1].real(), expected, 1e-8);
}

TEST(Tracker, TracksDecoupled2x2PathFromStartToTarget) {
  // G = { x^2 - 4, y^2 - 9 } (known start solution (2,3) at t=0)
  // F = { x^2 - 9, y^2 - 16 } (target solution (3,4) at t=1)
  // Decoupled (diagonal) system: each variable follows its own smooth
  // univariate path, same shape as the proven 1-variable case above.
  Variable x("x"), y("y");
  auto g1 = pow(x, 2) - 4.0;
  auto g2 = pow(y, 2) - 9.0;
  auto f1 = pow(x, 2) - 9.0;
  auto f2 = pow(y, 2) - 16.0;
  System<double> start({g1, g2}, {x, y});
  System<double> target({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(start, target);
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor, 1e-12, 20, 0.05, 1e-8, 0.1, 200);

  auto result = tracker.track({cd(2.0, 0.0), cd(3.0, 0.0)}, 0.0, 1.0);

  EXPECT_TRUE(result.success);
  EXPECT_NEAR(result.solution[0].real(), 3.0, 1e-8);
  EXPECT_NEAR(result.solution[1].real(), 4.0, 1e-8);
}

TEST(Tracker, FailsWhenStepBudgetExhausted) {
  // Same path as above, but capped at 1 step so it cannot possibly reach t=1
  Variable x("x");
  System<double> start({pow(x, 2) - 4.0}, {x});
  System<double> target({pow(x, 2) - 9.0}, {x});
  StraightLineHomotopy<double> H(start, target);
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor, 1e-12, 20, 0.1, 1e-8, 0.1, /*max_steps=*/1);

  auto result = tracker.track({cd(2.0, 0.0)}, 0.0, 1.0);

  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.steps, 1);
  EXPECT_NEAR(result.t, 0.1, 1e-9);
}
