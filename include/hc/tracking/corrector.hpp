#pragma once
#include <complex>
#include <vector>

#include <Eigen/Dense>

#include <hc/homotopy/homotopy.hpp>

namespace hc
{

struct CorrectorResult
{
  bool converged;
  int iterations;
};

// Newton's method: x <- x - Hx^-1 * H, evaluated at fixed t.
template <typename Scalar>
class NewtonCorrector
{
public:
  NewtonCorrector(Scalar tol, int max_iters) : tol_(tol), max_iters_(max_iters)
  {
  }

  CorrectorResult correct(std::vector<std::complex<Scalar>>& x,
                          const Homotopy<Scalar>& H, Scalar t) const
  {
    using Complex = std::complex<Scalar>;
    int n = H.num_equations();

    std::vector<Complex> H_val;
    std::vector<std::vector<Complex>> Hx_val;
    std::vector<Complex> Ht_val;

    for (int iter = 0; iter < max_iters_; ++iter)
    {
      H.evaluate_with_jacobian(x, t, H_val, Hx_val, Ht_val);

      Eigen::Matrix<Complex, Eigen::Dynamic, 1> r(n);
      for (int i = 0; i < n; ++i)
        r(i) = H_val[i];

      if (r.norm() < tol_)
        return {true, iter};

      Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic> J(n, n);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          J(i, j) = Hx_val[i][j];

      Eigen::Matrix<Complex, Eigen::Dynamic, 1> delta = J.partialPivLu().solve(r);

      for (int i = 0; i < n; ++i)
        x[i] -= delta(i);
    }

    H.evaluate(x, t, H_val);
    Eigen::Matrix<Complex, Eigen::Dynamic, 1> r(n);
    for (int i = 0; i < n; ++i)
      r(i) = H_val[i];

    return {r.norm() < tol_, max_iters_};
  }

private:
  Scalar tol_;
  int max_iters_;
};

} // namespace hc
