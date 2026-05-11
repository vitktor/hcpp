#include <gtest/gtest.h>
#include <hc/core/system.hpp>

using namespace hc;

TEST(System, ExplicitVars) {
  Variable x("x"), y("y");
  Polynomial<double> p({1.0, 1.0}, {{1, 0}, {0, 1}}, {x, y});
  System<double> s({p}, {x, y});
  EXPECT_EQ(s.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(s.getPolynomials().size(), 1u);
}

TEST(System, InferredVarsSinglePoly) {
  // x^3 - 2*x*y + y^2 - 1
  Variable x("x"), y("y");
  Polynomial<double> p({1.0, -2.0, 1.0, -1.0}, {{3, 0}, {1, 1}, {0, 2}, {0, 0}}, {x, y});
  System<double> s({p});
  EXPECT_EQ(s.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(s.getPolynomials()[0].degree(), 3);
  EXPECT_EQ(s.getPolynomials()[0].degree(x), 3);
  EXPECT_EQ(s.getPolynomials()[0].degree(y), 2);
}

TEST(System, InferredVarsOverlapping) {
  // p1: x^2 + x*y - 1  (uses {x, y})
  // p2: y^2 + y*z + z  (uses {y, z})
  // merged vars: {x, y, z}
  Variable x("x"), y("y"), z("z");
  Polynomial<double> p1({1.0, 1.0, -1.0}, {{2, 0}, {1, 1}, {0, 0}}, {x, y});
  Polynomial<double> p2({1.0, 1.0, 1.0}, {{2, 0}, {1, 1}, {0, 1}}, {y, z});
  System<double> s({p1, p2});
  EXPECT_EQ(s.getVariables(), (std::vector<Variable>{x, y, z}));
  EXPECT_EQ(s.getPolynomials().size(), 2u);
  EXPECT_EQ(s.getPolynomials()[0].degree(), 2);
  EXPECT_EQ(s.getPolynomials()[1].degree(), 2);
}

TEST(System, InferredVarsDisjoint) {
  // p1: x^2 - 1  (uses {x})
  // p2: y^3 + y  (uses {y})
  // merged vars: {x, y}
  Variable x("x"), y("y");
  Polynomial<double> p1({1.0, -1.0}, {{2}, {0}}, {x});
  Polynomial<double> p2({1.0, 1.0}, {{3}, {1}}, {y});
  System<double> s({p1, p2});
  EXPECT_EQ(s.getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(s.getPolynomials()[0].degree(x), 2);
  EXPECT_EQ(s.getPolynomials()[1].degree(y), 3);
}
