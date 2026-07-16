#include <gtest/gtest.h>
#include <hc/core/calculus.hpp>

using namespace hc;

// ---- differentiate ----

TEST(Differentiate, Constant) {
  // d/dx (5) = 0; the zero polynomial has no terms, and since the constant
  // 5 has no ambient variables to begin with, none are preserved either.
  Variable x("x");
  Polynomial<double> p(5.0);
  auto result = differentiate(p, x);
  EXPECT_EQ(result.degree(), -1);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{}));
  EXPECT_EQ(result.getExponents(), (std::vector<std::vector<int>>{}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{}));
}

TEST(Differentiate, SingleVariable) {
  // d/dx (x^3) = 3x^2
  Variable x("x");
  auto p = pow(x, 3);
  auto result = differentiate(p, x);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{3.0}));
  EXPECT_EQ(result.getExponents(), (std::vector<std::vector<int>>{{2}}));
  EXPECT_EQ(result.degree(x), 2);
}

TEST(Differentiate, OtherVariable) {
  // d/dy (x^3) = 0; the zero polynomial keeps x as its ambient variable
  Variable x("x"), y("y");
  auto p = pow(x, 3);
  auto result = differentiate(p, y);
  EXPECT_EQ(result.degree(), -1);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{}));
  EXPECT_EQ(result.getExponents(), (std::vector<std::vector<int>>{}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x}));
}

TEST(Differentiate, MultipleTerms) {
  // d/dx (x^3 + 2*x^2 + x + 1) = 3x^2 + 4x + 1
  Variable x("x");
  auto p = pow(x, 3) + 2.0*pow(x, 2) + x + 1.0;
  auto result = differentiate(p, x);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{3.0, 4.0, 1.0}));
  EXPECT_EQ(result.getExponents(), (std::vector<std::vector<int>>{{2}, {1}, {0}}));
  EXPECT_EQ(result.degree(x), 2);
}

TEST(Differentiate, MultipleVariables) {
  // d/dx (x^2 + x*y + y^2) = 2x + y
  Variable x("x"), y("y");
  auto p = pow(x, 2) + x*y + pow(y, 2);
  auto dx = differentiate(p, x);
  EXPECT_EQ(dx.getCoefficients(), (std::vector<double>{2.0, 1.0}));
  EXPECT_EQ(dx.getExponents(), (std::vector<std::vector<int>>{{1, 0}, {0, 1}}));
  EXPECT_EQ(dx.degree(x), 1);
  EXPECT_EQ(dx.degree(y), 1);
}

TEST(Differentiate, MultipleVariablesWrtY) {
  // d/dy (x^2 + x*y + y^2) = x + 2y
  Variable x("x"), y("y");
  auto p = pow(x, 2) + x*y + pow(y, 2);
  auto dy = differentiate(p, y);
  EXPECT_EQ(dy.getCoefficients(), (std::vector<double>{1.0, 2.0}));
  EXPECT_EQ(dy.getExponents(), (std::vector<std::vector<int>>{{1, 0}, {0, 1}}));
  EXPECT_EQ(dy.degree(x), 1);
  EXPECT_EQ(dy.degree(y), 1);
}

// ---- evaluate ----

TEST(Evaluate, Constant) {
  // p = 5, evaluate at x=2 -> 5
  Variable x("x");
  Polynomial<double> p(5.0);
  EXPECT_DOUBLE_EQ(evaluate(p, std::vector<double>{}), 5.0);
}

TEST(Evaluate, Linear) {
  // p = x + 1, evaluate at x=3 -> 4
  Variable x("x");
  auto p = x + 1.0;
  EXPECT_DOUBLE_EQ(evaluate(p, std::vector<double>{3.0}), 4.0);
}

TEST(Evaluate, Quadratic) {
  // p = x^2, evaluate at x=3 -> 9
  Variable x("x");
  auto p = pow(x, 2);
  EXPECT_DOUBLE_EQ(evaluate(p, std::vector<double>{3.0}), 9.0);
}

TEST(Evaluate, MultiVariable) {
  // p = x^2 + x*y + y^2, evaluate at (1, 2) -> 1 + 2 + 4 = 7
  Variable x("x"), y("y");
  auto p = pow(x, 2) + x*y + pow(y, 2);
  EXPECT_DOUBLE_EQ(evaluate(p, std::vector<double>{1.0, 2.0}), 7.0);
}

TEST(Evaluate, NegativeCoefficients) {
  // p = x^2 - 1, evaluate at x=1 -> 0
  Variable x("x");
  auto p = pow(x, 2) - 1.0;
  EXPECT_DOUBLE_EQ(evaluate(p, std::vector<double>{1.0}), 0.0);
}

TEST(Evaluate, ZeroPolynomialFromDifferentiate) {
  // Jacobian entries that come out to zero must remain safely evaluable
  // against a point sized for the polynomial's ambient variables.
  Variable x("x"), y("y");
  auto p = (pow(x, 2) - 1.0).update_vars({x, y}); // x^2 - 1, with y as an ambient variable
  auto dy = differentiate(p, y);
  EXPECT_DOUBLE_EQ(evaluate(dy, std::vector<double>{3.0, 4.0}), 0.0);
}
