#include <gtest/gtest.h>
#include <hc/core/polynomial.hpp>

using namespace hc;

TEST(Polynomial, Constructor) {
  Eigen::MatrixXi exps = Eigen::MatrixXi::Identity(1, 1);
  Polynomial<double> p({2.0}, exps, {Variable("x")});
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{Variable("x")}));
}

TEST(Polynomial, Degree) {
  // x^2 + y^3
  Eigen::MatrixXi exps(2, 2);
  exps << 2, 0, 0, 3;
  Polynomial<double> p({1.0, 1.0}, exps, {Variable("x"), Variable("y")});
  EXPECT_EQ(p.degree(), 3);
}

TEST(Polynomial, DegreeOfVariable) {
  // x^2 + y^3
  Eigen::MatrixXi exps(2, 2);
  exps << 2, 0, 0, 3;
  Polynomial<double> p({1.0, 1.0}, exps, {Variable("x"), Variable("y")});
  EXPECT_EQ(p.degree(Variable("x")), 2);
  EXPECT_EQ(p.degree(Variable("y")), 3);
  EXPECT_EQ(p.degree(Variable("z")), 0);
}

TEST(VariableAddition, SameVariable) {
  // x + x = 2x
  Variable x("x");
  auto p = x + x;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(p.degree(x), 1);
}

TEST(VariableAddition, DifferentVariables) {
  // x + y
  Variable x("x"), y("y");
  auto p = x + y;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0, 1.0}));
  EXPECT_EQ(p.degree(x), 1);
  EXPECT_EQ(p.degree(y), 1);
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
}

TEST(VariableAddition, DifferentVariablesSorted) {
  // y + x should produce the same sorted variable order as x + y
  Variable x("x"), y("y");
  auto p = y + x;
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
}

TEST(VariableNegation, Simple) {
  // -x
  Variable x("x");
  auto p = -x;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{-1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x}));
  EXPECT_EQ(p.degree(x), 1);
  EXPECT_EQ(p.degree(), 1);
}

TEST(VariableSubtraction, SameVariable) {
  // x - x = 0
  Variable x("x");
  auto p = x - x;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{0.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{}));
}

TEST(VariableSubtraction, DifferentVariables) {
  // x - y
  Variable x("x"), y("y");
  auto p = x - y;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0, -1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(p.degree(x), 1);
  EXPECT_EQ(p.degree(y), 1);
}

TEST(VariableSubtraction, DifferentVariablesSorted) {
  // y - x: vars sorted as {x, y}, coeffs {-1, 1}
  Variable x("x"), y("y");
  auto p = y - x;
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{-1.0, 1.0}));
}

TEST(VariableMultiplication, SameVariable) {
  // x * x = x^2
  Variable x("x");
  auto p = x * x;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x}));
  EXPECT_EQ(p.degree(x), 2);
  EXPECT_EQ(p.degree(), 2);
}

TEST(VariableMultiplication, DifferentVariables) {
  // x * y = xy
  Variable x("x"), y("y");
  auto p = x * y;
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(p.degree(x), 1);
  EXPECT_EQ(p.degree(y), 1);
  EXPECT_EQ(p.degree(), 2);
}

TEST(VariableMultiplication, DifferentVariablesSorted) {
  // y * x: vars sorted as {x, y}
  Variable x("x"), y("y");
  auto p = y * x;
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x, y}));
}

TEST(VariablePow, ZeroExponent) {
  // x^0 = 1 (constant, no variables)
  Variable x("x");
  auto p = pow(x, 0);
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{}));
  EXPECT_EQ(p.degree(), 0);
}

TEST(VariablePow, OneExponent) {
  // x^1 = x
  Variable x("x");
  auto p = pow(x, 1);
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x}));
  EXPECT_EQ(p.degree(x), 1);
}

TEST(VariablePow, HigherExponent) {
  // x^5
  Variable x("x");
  auto p = pow(x, 5);
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{x}));
  EXPECT_EQ(p.degree(x), 5);
  EXPECT_EQ(p.degree(), 5);
}
