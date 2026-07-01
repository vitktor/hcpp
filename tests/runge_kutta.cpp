#include <gtest/gtest.h>
#include <hc/homotopy/straight_line.hpp>
#include <hc/tracking/corrector.hpp>
#include <hc/tracking/runge_kutta.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(EulerPredictor, LinearHomotopyExactStep) {
  // G = x - 2, F = x - 5; H(x,t) = x - (2 + 3t), root(t) = 2 + 3t is linear in t,
  // so dx/dt = 3 everywhere and a single Euler step from t=1 to t=0 is exact.
  Variable x("x");
  System<double> start({Polynomial<double>({1.0, -2.0}, {{1}, {0}}, {x})}, {x});
  System<double> target({Polynomial<double>({1.0, -5.0}, {{1}, {0}}, {x})}, {x});
  StraightLineHomotopy<double> H(start, target);

  EulerPredictor<double> predictor;
  std::vector<cd> x0{cd(5.0, 0.0)}; // root at t=1
  std::vector<cd> x1;
  predictor.predict(x0, 1.0, -1.0, H, x1);

  EXPECT_NEAR(x1[0].real(), 2.0, 1e-12);
  EXPECT_NEAR(x1[0].imag(), 0.0, 1e-12);
}

TEST(EulerPredictor, PredictCorrectTracksNonlinearPath) {
  // G = x^2 - 4 (root 2), F = x^2 - 9 (root 3)
  // H(x,t) = x^2 - (4 + 5t); root(t) = sqrt(4+5t) is nonlinear in t, so a single
  // large Euler step would overshoot -- track in small predict+correct steps instead.
  Variable x("x");
  System<double> start({Polynomial<double>({1.0, -4.0}, {{2}, {0}}, {x})}, {x});
  System<double> target({Polynomial<double>({1.0, -9.0}, {{2}, {0}}, {x})}, {x});
  StraightLineHomotopy<double> H(start, target);

  EulerPredictor<double> predictor;
  NewtonCorrector<double> corrector(1e-12, 20);

  std::vector<cd> point{cd(3.0, 0.0)}; // root at t=1
  double t = 1.0;
  const double dt = -0.05;

  for (int i = 0; i < 20; ++i) {
    std::vector<cd> predicted;
    predictor.predict(point, t, dt, H, predicted);
    t += dt;
    auto result = corrector.correct(predicted, H, t);
    ASSERT_TRUE(result.converged);
    point = predicted;
  }

  EXPECT_NEAR(t, 0.0, 1e-9);
  EXPECT_NEAR(point[0].real(), 2.0, 1e-9);
  EXPECT_NEAR(point[0].imag(), 0.0, 1e-9);
}
