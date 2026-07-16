#include <cmath>

#include <gtest/gtest.h>
#include <hc/homotopy/straight_line.hpp>
#include <hc/tracking/corrector.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(NewtonCorrector, UnivariatePositiveRoot) {
  // F = x^2 - 1; start == target so H(x,t) = F(x) for every t
  Variable x("x");
  auto f = pow(x, 2) - 1.0;
  System<double> sys({f}, {x});
  StraightLineHomotopy<double> H(sys, sys);

  NewtonCorrector<double> corrector(1e-12, 20);
  std::vector<cd> point{cd(1.5, 0.0)};
  auto result = corrector.correct(point, H, 0.0);

  EXPECT_TRUE(result.converged);
  EXPECT_NEAR(point[0].real(), 1.0, 1e-9);
  EXPECT_NEAR(point[0].imag(), 0.0, 1e-9);
}

TEST(NewtonCorrector, UnivariateNegativeRoot) {
  // Starting closer to -1 should converge to the other root
  Variable x("x");
  auto f = pow(x, 2) - 1.0;
  System<double> sys({f}, {x});
  StraightLineHomotopy<double> H(sys, sys);

  NewtonCorrector<double> corrector(1e-12, 20);
  std::vector<cd> point{cd(-1.5, 0.0)};
  auto result = corrector.correct(point, H, 0.0);

  EXPECT_TRUE(result.converged);
  EXPECT_NEAR(point[0].real(), -1.0, 1e-9);
}

TEST(NewtonCorrector, CircleAndLine) {
  // F = { x^2 + y^2 - 1, x - y }; solution x=y=1/sqrt(2)
  Variable x("x"), y("y");
  auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
  auto f2 = x - y;
  System<double> sys({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(sys, sys);

  NewtonCorrector<double> corrector(1e-12, 20);
  std::vector<cd> point{cd(1.0, 0.0), cd(1.0, 0.0)};
  auto result = corrector.correct(point, H, 0.0);

  double expected = 1.0 / std::sqrt(2.0);
  EXPECT_TRUE(result.converged);
  EXPECT_NEAR(point[0].real(), expected, 1e-9);
  EXPECT_NEAR(point[1].real(), expected, 1e-9);
}

TEST(NewtonCorrector, AlreadyConverged) {
  // Starting exactly at the root should converge in 0 iterations
  Variable x("x");
  auto f = pow(x, 2) - 1.0;
  System<double> sys({f}, {x});
  StraightLineHomotopy<double> H(sys, sys);

  NewtonCorrector<double> corrector(1e-12, 20);
  std::vector<cd> point{cd(1.0, 0.0)};
  auto result = corrector.correct(point, H, 0.0);

  EXPECT_TRUE(result.converged);
  EXPECT_EQ(result.iterations, 0);
}
