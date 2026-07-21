#pragma once
#include <cmath>
#include <complex>
#include <iomanip>
#include <numbers>
#include <random>
#include <sstream>
#include <string>

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

// Formats a complex number as "a+bj" (or "a-bj"), unlike operator<<'s
// default "(a,b)".
template <typename Scalar>
std::string to_string(std::complex<Scalar> z, int precision = 6)
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << z.real()
      << (z.imag() < 0 ? "-" : "+") << std::abs(z.imag()) << "j";
  return oss.str();
}

} // namespace hc
