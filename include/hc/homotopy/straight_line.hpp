#pragma once
#include <complex>
#include <utility>
#include <vector>

#include <hc/core/calculus.hpp>
#include <hc/core/system.hpp>
#include <hc/homotopy/homotopy.hpp>
#include <hc/util.hpp>

namespace hc
{

// gamma defaults to a random unit-circle value (the "gamma trick", see
// Bates, Hauenstein, Sommese & Wampler, "Numerically Solving Polynomial
// Systems with Bertini"); pass gamma=1 explicitly for the plain homotopy.
template <typename Scalar>
struct StraightLineHomotopyOptions
{
  std::complex<Scalar> gamma = random_gamma<Scalar>();
};

// H(x,t) = gamma*(1-t)*start(x) + t*target(x). start and target must share
// the same variable list in the same order, since Hx_out columns are
// indexed by that shared order.
template <typename Scalar>
class StraightLineHomotopy : public Homotopy<Scalar>
{
public:
  StraightLineHomotopy(System<Scalar> start, System<Scalar> target,
                       StraightLineHomotopyOptions<Scalar> options = {})
      : start_(std::move(start)), target_(std::move(target)), gamma_(options.gamma)
  {
  }

  int num_equations() const override
  {
    return static_cast<int>(target_.getPolynomials().size());
  }

  int num_variables() const override
  {
    return static_cast<int>(target_.getVariables().size());
  }

  std::complex<Scalar> gamma() const { return gamma_; }

  void evaluate(const std::vector<std::complex<Scalar>>& x, Scalar t,
               std::vector<std::complex<Scalar>>& H_out) const override
  {
    const auto& G = start_.getPolynomials();
    const auto& F = target_.getPolynomials();
    std::complex<Scalar> gamma_one_minus_t = gamma_ * (Scalar(1) - t);
    std::complex<Scalar> tc(t);

    H_out.resize(F.size());
    for (size_t i = 0; i < F.size(); ++i)
      H_out[i] = gamma_one_minus_t * hc::evaluate(G[i], x) + tc * hc::evaluate(F[i], x);
  }

  void evaluate_with_jacobian(
      const std::vector<std::complex<Scalar>>& x, Scalar t,
      std::vector<std::complex<Scalar>>& H_out,
      std::vector<std::vector<std::complex<Scalar>>>& Hx_out,
      std::vector<std::complex<Scalar>>& Ht_out) const override
  {
    const auto& G = start_.getPolynomials();
    const auto& F = target_.getPolynomials();
    const auto& Gx = start_.getJacobian();
    const auto& Fx = target_.getJacobian();
    size_t n_eq = F.size();
    size_t n_vars = target_.getVariables().size();

    std::complex<Scalar> gamma_one_minus_t = gamma_ * (Scalar(1) - t);
    std::complex<Scalar> tc(t);

    H_out.resize(n_eq);
    Ht_out.resize(n_eq);
    Hx_out.assign(n_eq, std::vector<std::complex<Scalar>>(n_vars));

    for (size_t i = 0; i < n_eq; ++i)
    {
      auto g_val = hc::evaluate(G[i], x);
      auto f_val = hc::evaluate(F[i], x);
      H_out[i] = gamma_one_minus_t * g_val + tc * f_val;
      Ht_out[i] = f_val - gamma_ * g_val;
      for (size_t j = 0; j < n_vars; ++j)
        Hx_out[i][j] = gamma_one_minus_t * hc::evaluate(Gx[i][j], x) + tc * hc::evaluate(Fx[i][j], x);
    }
  }

private:
  System<Scalar> start_;
  System<Scalar> target_;
  std::complex<Scalar> gamma_;
};

} // namespace hc
