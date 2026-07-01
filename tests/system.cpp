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

TEST(Jacobian, CircleAndLine) {
  // F = { x^2 + y^2 - 1, x - y }
  // J = | 2x   2y |
  //     | 1    -1  |
  Variable x("x"), y("y");
  Polynomial<double> f1({1.0, 1.0, -1.0}, {{2, 0}, {0, 2}, {0, 0}}, {x, y});
  Polynomial<double> f2({1.0, -1.0}, {{1, 0}, {0, 1}}, {x, y});
  System<double> sys({f1, f2}, {x, y});
  const auto& J = sys.getJacobian();
  EXPECT_EQ(J.size(), 2u);
  EXPECT_EQ(J[0].size(), 2u);
  // J[0][0] = d/dx(x^2 + y^2 - 1) = 2x
  EXPECT_EQ(J[0][0].getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(J[0][0].getExponents(), (std::vector<std::vector<int>>{{1, 0}}));
  EXPECT_EQ(J[0][0].degree(x), 1);
  // J[0][1] = d/dy(x^2 + y^2 - 1) = 2y
  EXPECT_EQ(J[0][1].getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(J[0][1].getExponents(), (std::vector<std::vector<int>>{{0, 1}}));
  EXPECT_EQ(J[0][1].degree(y), 1);
  // J[1][0] = d/dx(x - y) = 1
  EXPECT_EQ(J[1][0].getCoefficients(), (std::vector<double>{1.0}));
  EXPECT_EQ(J[1][0].getExponents(), (std::vector<std::vector<int>>{{0, 0}}));
  EXPECT_EQ(J[1][0].degree(), 0);
  // J[1][1] = d/dy(x - y) = -1
  EXPECT_EQ(J[1][1].getCoefficients(), (std::vector<double>{-1.0}));
  EXPECT_EQ(J[1][1].getExponents(), (std::vector<std::vector<int>>{{0, 0}}));
  EXPECT_EQ(J[1][1].degree(), 0);
}

TEST(Jacobian, TotalDegree2x2) {
  // F = { x^2 - 1, y^2 - 1 }
  // J = | 2x  0  |
  //     | 0   2y |
  Variable x("x"), y("y");
  Polynomial<double> f1({1.0, -1.0}, {{2, 0}, {0, 0}}, {x, y});
  Polynomial<double> f2({1.0, -1.0}, {{0, 2}, {0, 0}}, {x, y});
  System<double> sys({f1, f2}, {x, y});
  const auto& J = sys.getJacobian();
  EXPECT_EQ(J[0][0].getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(J[0][0].getExponents(), (std::vector<std::vector<int>>{{1, 0}}));
  EXPECT_EQ(J[0][0].degree(x), 1);
  // J[0][1], J[1][0] are the zero polynomial: no terms, but the ambient
  // variables {x, y} are preserved so they remain safely evaluable.
  EXPECT_EQ(J[0][1].getCoefficients(), (std::vector<double>{}));
  EXPECT_EQ(J[0][1].getExponents(), (std::vector<std::vector<int>>{}));
  EXPECT_EQ(J[0][1].getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(J[0][1].degree(), -1);
  EXPECT_EQ(J[1][0].getCoefficients(), (std::vector<double>{}));
  EXPECT_EQ(J[1][0].getExponents(), (std::vector<std::vector<int>>{}));
  EXPECT_EQ(J[1][0].getVariables(), (std::vector<Variable>{x, y}));
  EXPECT_EQ(J[1][0].degree(), -1);
  EXPECT_EQ(J[1][1].getCoefficients(), (std::vector<double>{2.0}));
  EXPECT_EQ(J[1][1].getExponents(), (std::vector<std::vector<int>>{{0, 1}}));
  EXPECT_EQ(J[1][1].degree(y), 1);
}
