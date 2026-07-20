#include <gtest/gtest.h>
#include <hc/homotopy/straight_line.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(StraightLineHomotopy, Univariate) {
  // G = x - 2 (root 2), F = x - 5 (root 5)
  // H(x,t) = x - (2 + 3t)
  Variable x("x");
  System<double> start({x - 2.0}, {x});
  System<double> target({x - 5.0}, {x});
  StraightLineHomotopy<double> H(start, target);

  std::vector<cd> point{cd(2.0, 0.0)};
  std::vector<cd> H_out;

  H.evaluate(point, 0.0, H_out);
  EXPECT_DOUBLE_EQ(H_out[0].real(), 0.0);

  H.evaluate(point, 1.0, H_out);
  EXPECT_DOUBLE_EQ(H_out[0].real(), -3.0);
}

TEST(StraightLineHomotopy, EvaluateWithJacobianAtStart) {
  // G = {x^2 - 1, y^2 - 1}, F = {x^2 + y^2 - 1, x - y}
  Variable x("x"), y("y");
  auto g1 = pow(x, 2) - 1.0;
  auto g2 = pow(y, 2) - 1.0;
  auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
  auto f2 = x - y;
  System<double> start({g1, g2}, {x, y});
  System<double> target({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(start, target, {.gamma = cd(1.0, 0.0)});

  std::vector<cd> point{cd(1.0, 0.0), cd(1.0, 0.0)};
  std::vector<cd> H_out, Ht_out;
  std::vector<std::vector<cd>> Hx_out;

  H.evaluate_with_jacobian(point, 0.0, H_out, Hx_out, Ht_out);

  // H(x,0) = G(x): {1-1, 1-1} = {0, 0}
  EXPECT_DOUBLE_EQ(H_out[0].real(), 0.0);
  EXPECT_DOUBLE_EQ(H_out[1].real(), 0.0);

  // Hx(x,0) = Gx(x) = [[2x, 0], [0, 2y]] = [[2, 0], [0, 2]]
  EXPECT_DOUBLE_EQ(Hx_out[0][0].real(), 2.0);
  EXPECT_DOUBLE_EQ(Hx_out[0][1].real(), 0.0);
  EXPECT_DOUBLE_EQ(Hx_out[1][0].real(), 0.0);
  EXPECT_DOUBLE_EQ(Hx_out[1][1].real(), 2.0);

  // Ht(x) = F(x) - G(x) = {1-0, 0-0} = {1, 0}
  EXPECT_DOUBLE_EQ(Ht_out[0].real(), 1.0);
  EXPECT_DOUBLE_EQ(Ht_out[1].real(), 0.0);
}

TEST(StraightLineHomotopy, GammaDefaultsToUnitCircle) {
  Variable x("x");
  System<double> start({x - 2.0}, {x});
  System<double> target({x - 5.0}, {x});
  StraightLineHomotopy<double> H(start, target);
  EXPECT_NEAR(std::abs(H.gamma()), 1.0, 1e-12);
}

TEST(StraightLineHomotopy, GammaTrick) {
  // G = x - 2, F = x - 5, gamma = i
  // H(x,t) = i*(1-t)*(x-2) + t*(x-5)
  Variable x("x");
  System<double> start({x - 2.0}, {x});
  System<double> target({x - 5.0}, {x});
  cd gamma(0.0, 1.0);
  StraightLineHomotopy<double> H(start, target, {.gamma = gamma});
  EXPECT_EQ(H.gamma(), gamma);

  std::vector<cd> point{cd(0.0, 0.0)};
  std::vector<cd> H_out, Ht_out;
  std::vector<std::vector<cd>> Hx_out;

  // At t=0: H = gamma*G(x) = i*(0-2) = -2i
  H.evaluate_with_jacobian(point, 0.0, H_out, Hx_out, Ht_out);
  EXPECT_NEAR(H_out[0].real(), 0.0, 1e-12);
  EXPECT_NEAR(H_out[0].imag(), -2.0, 1e-12);
  // Hx(x,0) = gamma*Gx(x) = i*1 = i
  EXPECT_NEAR(Hx_out[0][0].real(), 0.0, 1e-12);
  EXPECT_NEAR(Hx_out[0][0].imag(), 1.0, 1e-12);
  // Ht(x) = F(x) - gamma*G(x) = (0-5) - i*(0-2) = -5 + 2i
  EXPECT_NEAR(Ht_out[0].real(), -5.0, 1e-12);
  EXPECT_NEAR(Ht_out[0].imag(), 2.0, 1e-12);

  // At t=1, H = F(x) regardless of gamma
  H.evaluate(point, 1.0, H_out);
  EXPECT_NEAR(H_out[0].real(), -5.0, 1e-12);
  EXPECT_NEAR(H_out[0].imag(), 0.0, 1e-12);
}

TEST(StraightLineHomotopy, EvaluateWithJacobianAtTarget) {
  Variable x("x"), y("y");
  auto g1 = pow(x, 2) - 1.0;
  auto g2 = pow(y, 2) - 1.0;
  auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
  auto f2 = x - y;
  System<double> start({g1, g2}, {x, y});
  System<double> target({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(start, target);

  std::vector<cd> point{cd(1.0, 0.0), cd(1.0, 0.0)};
  std::vector<cd> H_out, Ht_out;
  std::vector<std::vector<cd>> Hx_out;

  H.evaluate_with_jacobian(point, 1.0, H_out, Hx_out, Ht_out);

  // H(x,1) = F(x): {1+1-1, 1-1} = {1, 0}
  EXPECT_DOUBLE_EQ(H_out[0].real(), 1.0);
  EXPECT_DOUBLE_EQ(H_out[1].real(), 0.0);

  // Hx(x,1) = Fx(x) = [[2x, 2y], [1, -1]] = [[2, 2], [1, -1]]
  EXPECT_DOUBLE_EQ(Hx_out[0][0].real(), 2.0);
  EXPECT_DOUBLE_EQ(Hx_out[0][1].real(), 2.0);
  EXPECT_DOUBLE_EQ(Hx_out[1][0].real(), 1.0);
  EXPECT_DOUBLE_EQ(Hx_out[1][1].real(), -1.0);
}
