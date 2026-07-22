#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>
#include <hc/solver/solver.hpp>

using namespace hc;

TEST(Solver, UnivariateQuadratic) {
  // x^2 - 9 = 0, roots +-3. Bezout number (2) matches the real root count
  // exactly, so both paths should converge.
  Variable x("x");
  System<double> target({pow(x, 2) - 9.0}, {x});

  Solver<double> solver;
  auto result = solver.solve(target);

  auto sols = result.real_solutions();
  ASSERT_EQ(sols.size(), 2u);
  std::vector<double> real_parts;
  for (const auto& s : sols)
    real_parts.push_back(s[0].real());
  std::sort(real_parts.begin(), real_parts.end());
  EXPECT_NEAR(real_parts[0], -3.0, 1e-8);
  EXPECT_NEAR(real_parts[1], 3.0, 1e-8);
}

TEST(Solver, TwoVariableDecoupled) {
  // { x^2 - 9, y^2 - 16 } = 0, roots (+-3, +-4). Bezout number 4 matches
  // the real root count exactly.
  Variable x("x"), y("y");
  System<double> target({pow(x, 2) - 9.0, pow(y, 2) - 16.0}, {x, y});

  Solver<double> solver;
  auto result = solver.solve(target);

  auto sols = result.real_solutions();
  EXPECT_EQ(sols.size(), 4u);
  for (const auto& s : sols) {
    EXPECT_NEAR(std::abs(s[0].real()), 3.0, 1e-7);
    EXPECT_NEAR(std::abs(s[1].real()), 4.0, 1e-7);
  }
}

TEST(Solver, OptionsOverrideStepSize) {
  Variable x("x");
  System<double> target({pow(x, 2) - 9.0}, {x});

  Solver<double> solver({.dt_max = 0.2, .max_steps = 50});
  auto result = solver.solve(target);

  EXPECT_EQ(result.real_solutions().size(), 2u);
}
