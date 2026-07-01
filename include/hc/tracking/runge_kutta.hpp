#pragma once
#include <complex>
#include <vector>

#include <Eigen/Dense>

#include <hc/homotopy/homotopy.hpp>
#include <hc/tracking/predictor.hpp>

namespace hc
{

// Euler step: x <- x + dt * dx/dt, where dx/dt = -Hx^-1 * Ht at (x, t).
template <typename Scalar>
class EulerPredictor : public Predictor<Scalar>
{
public:
  void predict(const std::vector<std::complex<Scalar>>& x, Scalar t, Scalar dt,
              const Homotopy<Scalar>& H,
              std::vector<std::complex<Scalar>>& x_out) const override
  {
    using Complex = std::complex<Scalar>;
    int n = H.num_equations();

    std::vector<Complex> H_val;
    std::vector<std::vector<Complex>> Hx_val;
    std::vector<Complex> Ht_val;
    H.evaluate_with_jacobian(x, t, H_val, Hx_val, Ht_val);

    Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> J(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        J(i, j) = Hx_val[i][j];

    Eigen::Matrix<Complex, Eigen::Dynamic, 1> ht(n);
    for (int i = 0; i < n; ++i)
      ht(i) = -Ht_val[i];

    Eigen::Matrix<Complex, Eigen::Dynamic, 1> dxdt = J.partialPivLu().solve(ht);

    x_out.resize(n);
    for (int i = 0; i < n; ++i)
      x_out[i] = x[i] + Complex(dt) * dxdt(i);
  }
};

} // namespace hc
