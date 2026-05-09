#include <gtest/gtest.h>
#include <hc/core/polynomial.hpp>

using namespace hc;

TEST(Polynomial, Constructor) {
  Polynomial<double> p({2.0}, {{1}}, {Variable("x")});
  EXPECT_EQ(p.getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(p.getVariables(), (std::vector<Variable>{Variable("x")}));
}

TEST(Polynomial, Degree) {
  // x^2 + y^3: monomial 0 = x^2 (exps {2,0}), monomial 1 = y^3 (exps {0,3})
  Polynomial<double> p({1.0, 1.0}, {{2, 0}, {0, 3}}, {Variable("x"), Variable("y")});
  EXPECT_EQ(p.degree(), 3);
}

TEST(Polynomial, DegreeOfVariable) {
  // x^2 + y^3
  Polynomial<double> p({1.0, 1.0}, {{2, 0}, {0, 3}}, {Variable("x"), Variable("y")});
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

TEST(PolynomialNegation, Simple) {
  // -(x^2 + 2y) = -x^2 - 2y; decreasing lex puts x^2 before y
  Variable x("x"), y("y");
  Polynomial<double> p({1.0, 2.0}, {{2, 0}, {0, 1}}, {x, y});
  auto result = -p;
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{-1.0, -2.0}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x, y}));
}

TEST(PolynomialAddition, SameTerms) {
  // 2x + 3x = 5x
  Variable x("x");
  Polynomial<double> p1({2.0}, {{1}}, {x});
  Polynomial<double> p2({3.0}, {{1}}, {x});
  auto result = p1 + p2;
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{5.0}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x}));
}

TEST(PolynomialAddition, CancellingTerms) {
  // x + (-x) = 0
  Variable x("x");
  Polynomial<double> p1({1.0}, {{1}}, {x});
  Polynomial<double> p2({-1.0}, {{1}}, {x});
  auto result = p1 + p2;
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{}));
}

TEST(PolynomialAddition, DifferentVariables) {
  // x + y: degree 1 in both
  Variable x("x"), y("y");
  Polynomial<double> p1({1.0}, {{1}}, {x});
  Polynomial<double> p2({1.0}, {{1}}, {y});
  auto result = p1 + p2;
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(result.degree(x), 1);
  EXPECT_EQ(result.degree(y), 1);
}

TEST(PolynomialAddition, MergedVariables) {
  // (x + y) + (y + z) = x + 2y + z
  Variable x("x"), y("y"), z("z");
  Polynomial<double> p1({1.0, 1.0}, {{1, 0}, {0, 1}}, {x, y});
  Polynomial<double> p2({1.0, 1.0}, {{1, 0}, {0, 1}}, {y, z});
  auto result = p1 + p2;
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x, y, z}));
  EXPECT_EQ(result.degree(x), 1);
  EXPECT_EQ(result.degree(y), 1);
  EXPECT_EQ(result.degree(z), 1);
}

TEST(PolynomialMultiplication, SameVariable) {
  // x * x = x^2
  Variable x("x");
  Polynomial<double> p1({1.0}, {{1}}, {x});
  Polynomial<double> p2({1.0}, {{1}}, {x});
  auto result = p1 * p2;
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x}));
  EXPECT_EQ(result.degree(x), 2);
}

TEST(PolynomialMultiplication, DifferentVariables) {
  // (2x) * (3y) = 6xy
  Variable x("x"), y("y");
  Polynomial<double> p1({2.0}, {{1, 0}}, {x, y});
  Polynomial<double> p2({3.0}, {{0, 1}}, {x, y});
  auto result = p1 * p2;
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{6.0}));
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(result.degree(x), 1);
  EXPECT_EQ(result.degree(y), 1);
}

TEST(PolynomialMultiplication, LikeTerms) {
  // (x + y) * (x - y) = x^2 - y^2
  Variable x("x"), y("y");
  Polynomial<double> p1({1.0, 1.0}, {{1, 0}, {0, 1}}, {x, y});
  Polynomial<double> p2({1.0, -1.0}, {{1, 0}, {0, 1}}, {x, y});
  auto result = p1 * p2;
  EXPECT_EQ(result.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(result.degree(x), 2);
  EXPECT_EQ(result.degree(y), 2);
  // x^2 before y^2 in decreasing lex
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{1.0, -1.0}));
}

TEST(PolynomialMultiplication, CancellingTerms) {
  // (x + y) * (x - y) middle terms cancel; x^2 - xy + xy - y^2 = x^2 - y^2
  // already tested above; verify no xy term
  Variable x("x"), y("y");
  Polynomial<double> p1({1.0, 1.0}, {{1, 0}, {0, 1}}, {x, y});
  Polynomial<double> p2({1.0, -1.0}, {{1, 0}, {0, 1}}, {x, y});
  auto result = p1 * p2;
  EXPECT_EQ(result.getCoefficients().size(), 2u);
}
