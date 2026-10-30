#include <gtest/gtest.h>
#include <hc/solver/monodromy.hpp>

using namespace hc;
using cd = std::complex<double>;

namespace {

double residual(const System<double>& F, const std::vector<Variable>& params,
                const std::vector<cd>& x, const std::vector<cd>& p) {
  ParameterHomotopy<double> H(F, params, p, p);
  std::vector<cd> H_out;
  H.evaluate(x, 0.0, H_out);
  double r = 0.0;
  for (const auto& h : H_out)
    r = std::max(r, std::abs(h));
  return r;
}

void expect_distinct(const std::vector<std::vector<cd>>& sols) {
  for (size_t i = 0; i < sols.size(); ++i)
    for (size_t j = i + 1; j < sols.size(); ++j) {
      double dist = 0.0;
      for (size_t k = 0; k < sols[i].size(); ++k)
        dist = std::max(dist, std::abs(sols[i][k] - sols[j][k]));
      EXPECT_GT(dist, 1e-4);
    }
}

} // namespace

TEST(Monodromy, QuarticTrinomial) {
  // x^4 + a*x + b has 4 roots for generic (a, b). Start pair: pick x0 and a,
  // then b = -(x0^4 + a*x0) makes x0 a root.
  Variable x("x"), a("a"), b("b");
  System<double> F({pow(x, 4) + a * x + b});
  std::vector<Variable> params{a, b};

  cd x0(0.3, 0.7), a0(-0.5, 1.1);
  cd b0 = -(std::pow(x0, 4) + a0 * x0);

  MonodromySolver<double> solver;
  auto result = solver.solve(F, params, {x0}, {a0, b0});

  ASSERT_EQ(result.solutions.size(), 4u);
  for (const auto& s : result.solutions)
    EXPECT_LT(residual(F, params, s, result.parameters), 1e-10);
  expect_distinct(result.solutions);
}

TEST(Monodromy, TwoVariablesWithTargetCount) {
  // { x^2 + y^2 - 1 - a, x*y - b } has 4 solutions for generic (a, b).
  Variable x("x"), y("y"), a("a"), b("b");
  System<double> F({pow(x, 2) + pow(y, 2) - 1.0 - a, x * y - b});
  std::vector<Variable> params{a, b};

  cd x0(0.4, -0.2), y0(-0.9, 0.6);
  cd a0 = x0 * x0 + y0 * y0 - 1.0;
  cd b0 = x0 * y0;

  MonodromySolver<double> solver({.target_solutions_count = 4});
  auto result = solver.solve(F, params, {x0, y0}, {a0, b0});

  ASSERT_EQ(result.solutions.size(), 4u);
  for (const auto& s : result.solutions)
    EXPECT_LT(residual(F, params, s, result.parameters), 1e-10);
  expect_distinct(result.solutions);
}

TEST(Monodromy, StopsAfterLoopsWithoutProgress) {
  // x - a has a single solution, so no loop can ever find a new one.
  Variable x("x"), a("a");
  System<double> F({x - a});

  MonodromySolver<double> solver({.max_loops_no_progress = 3});
  auto result = solver.solve(F, {a}, {cd(0.5, 0.5)}, {cd(0.5, 0.5)});

  EXPECT_EQ(result.solutions.size(), 1u);
  EXPECT_EQ(result.loops, 3);
}
