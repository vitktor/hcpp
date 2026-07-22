#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>
#include <hc/homotopy/straight_line.hpp>
#include <hc/start_systems/total_degree.hpp>
#include <hc/tracking/runge_kutta.hpp>
#include <hc/tracking/tracker.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(TotalDegreeSystem, UnivariateQuadratic) {
  // target = x^2 - 9 (degree 2) -> G = x^2 - 1, start solutions {1, -1}
  Variable x("x");
  System<double> target({pow(x, 2) - 9.0}, {x});
  TotalDegreeSystem<double> start(target);

  const auto& g = start.system().getPolynomials();
  ASSERT_EQ(g.size(), 1u);
  EXPECT_EQ(g[0].getCoefficients(), (std::vector<double>{1.0, -1.0}));
  EXPECT_EQ(g[0].getExponents(), (std::vector<std::vector<int>>{{2}, {0}}));

  const auto& sols = start.start_solutions();
  ASSERT_EQ(sols.size(), 2u);
  EXPECT_NEAR(sols[0][0].real(), 1.0, 1e-12);
  EXPECT_NEAR(sols[0][0].imag(), 0.0, 1e-12);
  EXPECT_NEAR(sols[1][0].real(), -1.0, 1e-12);
  EXPECT_NEAR(sols[1][0].imag(), 0.0, 1e-12);
}

TEST(TotalDegreeSystem, TwoVariableMixedDegrees) {
  // target = {x^2 - 9, y^3 - 8}, degrees {2,3} -> G = {x^2 - 1, y^3 - 1}
  // Bezout number 2*3=6, Cartesian product of the 2nd and 3rd roots of unity.
  Variable x("x"), y("y");
  System<double> target({pow(x, 2) - 9.0, pow(y, 3) - 8.0}, {x, y});
  TotalDegreeSystem<double> start(target);

  const auto& g = start.system().getPolynomials();
  ASSERT_EQ(g.size(), 2u);
  EXPECT_EQ(g[0].getCoefficients(), (std::vector<double>{1.0, -1.0}));
  EXPECT_EQ(g[1].getCoefficients(), (std::vector<double>{1.0, -1.0}));

  const auto& sols = start.start_solutions();
  ASSERT_EQ(sols.size(), 6u);
  for (const auto& sol : sols) {
    ASSERT_EQ(sol.size(), 2u);
    // x-component is a 2nd root of unity: x^2 == 1
    EXPECT_NEAR(std::abs(sol[0] * sol[0] - cd(1.0, 0.0)), 0.0, 1e-12);
    // y-component is a 3rd root of unity: y^3 == 1
    EXPECT_NEAR(std::abs(sol[1] * sol[1] * sol[1] - cd(1.0, 0.0)), 0.0, 1e-12);
  }
}

TEST(TotalDegreeSystem, TracksToRealTarget) {
  // target = x^2 - 9 (roots +-3); Bezout number of G matches exactly (2),
  // so both start solutions should converge, one to each target root.
  Variable x("x");
  System<double> target({pow(x, 2) - 9.0}, {x});
  TotalDegreeSystem<double> start(target);

  StraightLineHomotopy<double> H(start.system(), target);
  EulerPredictor<double> predictor;
  Tracker<double> tracker(H, predictor, 1e-12, 20, 0.05, 1e-9, 0.1, 1000);

  std::vector<double> real_parts;
  for (const auto& x0 : start.start_solutions()) {
    auto result = tracker.track(x0, 0.0, 1.0);
    ASSERT_TRUE(result.success);
    EXPECT_NEAR(result.solution[0].imag(), 0.0, 1e-8);
    real_parts.push_back(result.solution[0].real());
  }
  std::sort(real_parts.begin(), real_parts.end());
  ASSERT_EQ(real_parts.size(), 2u);
  EXPECT_NEAR(real_parts[0], -3.0, 1e-8);
  EXPECT_NEAR(real_parts[1], 3.0, 1e-8);
}
