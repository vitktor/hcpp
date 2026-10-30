#include <gtest/gtest.h>
#include <hc/homotopy/parameter.hpp>
#include <hc/tracking/runge_kutta.hpp>
#include <hc/tracking/tracker.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(ParameterHomotopy, EvaluateWithJacobian) {
  // F = x^2 + a*x + b, parameters {a, b}. System vars are sorted (a, b, x),
  // so the parameters come before the unknown.
  Variable x("x"), a("a"), b("b");
  System<double> F({pow(x, 2) + a * x + b});
  ParameterHomotopy<double> H(F, {a, b}, {cd(1.0), cd(1.0)}, {cd(3.0), cd(5.0)});
  EXPECT_EQ(H.num_variables(), 1);
  EXPECT_EQ(H.num_parameters(), 2);

  std::vector<cd> point{cd(2.0)};
  std::vector<cd> H_out, Ht_out;
  std::vector<std::vector<cd>> Hx_out;

  // At t=0.5, p = (2, 3): H = 4 + 2*2 + 3 = 11
  H.evaluate_with_jacobian(point, 0.5, H_out, Hx_out, Ht_out);
  EXPECT_NEAR(H_out[0].real(), 11.0, 1e-12);
  // Hx = 2x + a = 4 + 2 = 6
  ASSERT_EQ(Hx_out[0].size(), 1u);
  EXPECT_NEAR(Hx_out[0][0].real(), 6.0, 1e-12);
  // Ht = dF/da * (3-1) + dF/db * (5-1) = 2*2 + 1*4 = 8
  EXPECT_NEAR(Ht_out[0].real(), 8.0, 1e-12);

  H.evaluate(point, 0.5, H_out);
  EXPECT_NEAR(H_out[0].real(), 11.0, 1e-12);
}

TEST(ParameterHomotopy, SetParameters) {
  Variable x("x"), a("a");
  System<double> F({pow(x, 2) - a});
  ParameterHomotopy<double> H(F, {a}, {cd(1.0)}, {cd(1.0)});
  H.set_parameters({cd(4.0)}, {cd(9.0)});

  std::vector<cd> H_out;
  H.evaluate({cd(3.0)}, 1.0, H_out);
  EXPECT_NEAR(std::abs(H_out[0]), 0.0, 1e-12);

  EXPECT_THROW(H.set_parameters({cd(1.0), cd(2.0)}, {cd(1.0)}), std::invalid_argument);
}

TEST(ParameterHomotopy, UnknownParameterThrows) {
  Variable x("x"), a("a"), c("c");
  System<double> F({pow(x, 2) - a});
  EXPECT_THROW(ParameterHomotopy<double>(F, {c}, {cd(1.0)}, {cd(1.0)}),
               std::invalid_argument);
}

TEST(ParameterHomotopy, TracksSolutionBetweenParameters) {
  // x^2 - a: track the root x=2 at a=4 to x=3 at a=9.
  Variable x("x"), a("a");
  System<double> F({pow(x, 2) - a});
  ParameterHomotopy<double> H(F, {a}, {cd(4.0)}, {cd(9.0)});
  EulerPredictor<double> predictor;
  Tracker<double> tracker(H, predictor, 1e-12, 20, 0.05, 1e-9, 0.1, 1000);

  auto result = tracker.track({cd(2.0)}, 0.0, 1.0);
  EXPECT_TRUE(result.success);
  EXPECT_NEAR(result.solution[0].real(), 3.0, 1e-8);
  EXPECT_NEAR(result.solution[0].imag(), 0.0, 1e-8);
}
