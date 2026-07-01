#pragma once
#include <complex>
#include <vector>

namespace hc
{

template <typename Scalar>
class Homotopy
{
public:
  virtual ~Homotopy() = default;

  virtual int num_equations() const = 0;
  virtual int num_variables() const = 0;

  virtual void evaluate(const std::vector<std::complex<Scalar>>& x, Scalar t,
                        std::vector<std::complex<Scalar>>& H_out) const = 0;

  virtual void evaluate_with_jacobian(
      const std::vector<std::complex<Scalar>>& x, Scalar t,
      std::vector<std::complex<Scalar>>& H_out,
      std::vector<std::vector<std::complex<Scalar>>>& Hx_out,
      std::vector<std::complex<Scalar>>& Ht_out) const = 0;
};

} // namespace hc
