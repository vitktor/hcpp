#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <random>

namespace hc
{

// Random point on the unit circle; fixed seed for reproducibility.
template <typename Scalar>
std::complex<Scalar> random_gamma()
{
  static std::mt19937 rng(12345);
  static std::uniform_real_distribution<Scalar> angle(Scalar(0),
                                                       Scalar(2) * std::numbers::pi_v<Scalar>);
  Scalar theta = angle(rng);
  return std::complex<Scalar>(std::cos(theta), std::sin(theta));
}

} // namespace hc
