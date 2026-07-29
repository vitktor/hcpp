#include <boost/multiprecision/cpp_bin_float.hpp>
#include <gtest/gtest.h>

#include <hc/core/calculus.hpp>
#include <hc/core/polynomial.hpp>

using namespace hc;
using HP = boost::multiprecision::cpp_bin_float_50;

// Demonstrates that hc::Polynomial/differentiate/evaluate -- the core,
// Scalar-templated layer of hcpp -- work unmodified with an arbitrary-
// precision Boost.Multiprecision type, not just double/long double. This
// is the reason Scalar is a template parameter throughout the codebase:
// re-tracking ill-conditioned paths at higher precision without rewriting
// algorithms.
//
// Stays on the real-valued core layer deliberately: std::complex<HP> hits
// a hard wall on this toolchain -- libc++'s std::complex<T>::abs() calls
// std::hypot() internally, which is gated to is_arithmetic types, so
// Eigen's partialPivLu() pivoting (used by NewtonCorrector/Tracker) can't
// take Scalar=HP without a dedicated multiprecision complex backend (e.g.
// Boost's mpc_complex, which needs GMP/MPFR/MPC) -- out of scope for a
// minimal first Boost dependency.
TEST(Multiprecision, NewtonsMethodBeyondDoublePrecision) {
  // F = x^2 - 2, root sqrt(2). Real-valued Newton's method built directly
  // from hc::evaluate/differentiate; built via the Polynomial constructor
  // directly since the `pow`/operator sugar (polynomial.hpp) is only
  // instantiated for double.
  Variable x("x");
  Polynomial<HP> f({HP(1), HP(-2)}, {{2}, {0}}, {x});
  Polynomial<HP> df = differentiate(f, x);

  HP guess("1.4");
  for (int i = 0; i < 12; ++i)
  {
    HP fx = evaluate(f, std::vector<HP>{guess});
    HP dfx = evaluate(df, std::vector<HP>{guess});
    guess -= fx / dfx;
  }

  // 1e-40 accuracy is far beyond double's ~1e-16 machine epsilon --
  // reachable only because Scalar=HP carries ~50 decimal digits here.
  HP sqrt2 = boost::multiprecision::sqrt(HP(2));
  EXPECT_LT(abs(guess - sqrt2), HP("1e-40"));
}
