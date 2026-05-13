#include <gtest/gtest.h>
#include <hc/core/differentiation.hpp>

using namespace hc;

TEST(Differentiate, Constant) {
  // d/dx (5) = 0
  Variable x("x");
  Polynomial<double> p(5.0);
  auto result = differentiate(p, x);
  EXPECT_EQ(result.degree(), 0);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{0.0}));
}

TEST(Differentiate, SingleVariable) {
  // d/dx (x^3) = 3x^2
  Variable x("x");
  Polynomial<double> p({1.0}, {{3}}, {x});
  auto result = differentiate(p, x);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{3.0}));
  EXPECT_EQ(result.degree(x), 2);
}

TEST(Differentiate, OtherVariable) {
  // d/dy (x^3) = 0
  Variable x("x"), y("y");
  Polynomial<double> p({1.0}, {{3}}, {x});
  auto result = differentiate(p, y);
  EXPECT_EQ(result.degree(), 0);
}

TEST(Differentiate, MultipleTerms) {
  // d/dx (x^3 + 2*x^2 + x + 1) = 3x^2 + 4x + 1
  Variable x("x");
  Polynomial<double> p({1.0, 2.0, 1.0, 1.0}, {{3}, {2}, {1}, {0}}, {x});
  auto result = differentiate(p, x);
  EXPECT_EQ(result.getCoefficients(), (std::vector<double>{3.0, 4.0, 1.0}));
  EXPECT_EQ(result.degree(x), 2);
}

TEST(Differentiate, MultipleVariables) {
  // d/dx (x^2 + x*y + y^2) = 2x + y
  Variable x("x"), y("y");
  Polynomial<double> p({1.0, 1.0, 1.0}, {{2, 0}, {1, 1}, {0, 2}}, {x, y});
  auto dx = differentiate(p, x);
  EXPECT_EQ(dx.getCoefficients(), (std::vector<double>{2.0, 1.0}));
  EXPECT_EQ(dx.degree(x), 1);
  EXPECT_EQ(dx.degree(y), 1);
}

TEST(Differentiate, MultipleVariablesWrtY) {
  // d/dy (x^2 + x*y + y^2) = x + 2y
  Variable x("x"), y("y");
  Polynomial<double> p({1.0, 1.0, 1.0}, {{2, 0}, {1, 1}, {0, 2}}, {x, y});
  auto dy = differentiate(p, y);
  EXPECT_EQ(dy.getCoefficients(), (std::vector<double>{1.0, 2.0}));
  EXPECT_EQ(dy.degree(x), 1);
  EXPECT_EQ(dy.degree(y), 1);
}
