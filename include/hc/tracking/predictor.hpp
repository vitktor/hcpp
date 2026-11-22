#pragma once
#include <complex>
#include <vector>

#include <hc/homotopy/homotopy.hpp>

namespace hc
{

// Thread safety: predict() is called concurrently on one shared instance,
// so implementations must not modify shared state (see Homotopy).
template <typename Scalar>
class Predictor
{
public:
  virtual ~Predictor() = default;

  // Predicts x at t + dt from the current point x at t.
  virtual void predict(const std::vector<std::complex<Scalar>>& x, Scalar t, Scalar dt,
                       const Homotopy<Scalar>& H,
                       std::vector<std::complex<Scalar>>& x_out) const = 0;
};

} // namespace hc
