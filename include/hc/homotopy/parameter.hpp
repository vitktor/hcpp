#pragma once
#include <algorithm>
#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include <hc/core/calculus.hpp>
#include <hc/core/system.hpp>
#include <hc/homotopy/homotopy.hpp>

namespace hc
{

// H(x,t) = F(x, p(t)), p(t) = (1-t)*p_start + t*p_target, for a family F
// given as one System over unknowns and parameters together. Because
// System keeps its variables sorted, unknowns and parameters may interleave;
// the unknowns x are F's variables minus `parameters`, in F's (sorted)
// order, and p_start/p_target are ordered like `parameters`.
//
// The parameter endpoints can be changed with set_parameters(), so one
// ParameterHomotopy (and a Tracker referencing it) can be reused for many
// segments, e.g. the edges of a monodromy loop.
template <typename Scalar>
class ParameterHomotopy : public Homotopy<Scalar>
{
public:
  using Complex = std::complex<Scalar>;

  ParameterHomotopy(System<Scalar> F, const std::vector<Variable>& parameters,
                    std::vector<Complex> p_start, std::vector<Complex> p_target)
      : F_(std::move(F))
  {
    const auto& vars = F_.getVariables();
    for (const auto& p : parameters)
    {
      auto it = std::find(vars.begin(), vars.end(), p);
      if (it == vars.end())
        throw std::invalid_argument("ParameterHomotopy: parameter " + p.name +
                                    " is not a variable of the system");
      param_idx_.push_back(static_cast<int>(it - vars.begin()));
    }
    for (size_t j = 0; j < vars.size(); ++j)
      if (std::find(parameters.begin(), parameters.end(), vars[j]) == parameters.end())
        var_idx_.push_back(static_cast<int>(j));

    set_parameters(std::move(p_start), std::move(p_target));
  }

  int num_equations() const override
  {
    return static_cast<int>(F_.getPolynomials().size());
  }

  int num_variables() const override { return static_cast<int>(var_idx_.size()); }

  int num_parameters() const { return static_cast<int>(param_idx_.size()); }

  void set_parameters(std::vector<Complex> p_start, std::vector<Complex> p_target)
  {
    if (p_start.size() != param_idx_.size() || p_target.size() != param_idx_.size())
      throw std::invalid_argument(
          "ParameterHomotopy: parameter values must match the number of parameters");
    p_start_ = std::move(p_start);
    p_target_ = std::move(p_target);
  }

  const std::vector<Complex>& start_parameters() const { return p_start_; }
  const std::vector<Complex>& target_parameters() const { return p_target_; }

  void evaluate(const std::vector<Complex>& x, Scalar t,
                std::vector<Complex>& H_out) const override
  {
    const auto& F = F_.getPolynomials();
    auto point = full_point(x, t);

    H_out.resize(F.size());
    for (size_t i = 0; i < F.size(); ++i)
      H_out[i] = hc::evaluate(F[i], point);
  }

  void evaluate_with_jacobian(const std::vector<Complex>& x, Scalar t,
                              std::vector<Complex>& H_out,
                              std::vector<std::vector<Complex>>& Hx_out,
                              std::vector<Complex>& Ht_out) const override
  {
    const auto& F = F_.getPolynomials();
    const auto& J = F_.getJacobian();
    auto point = full_point(x, t);
    size_t n_eq = F.size();

    H_out.resize(n_eq);
    Ht_out.resize(n_eq);
    Hx_out.assign(n_eq, std::vector<Complex>(var_idx_.size()));

    // dp/dt = p_target - p_start, so Ht = sum_k dF/dp_k * (p_target - p_start)_k.
    for (size_t i = 0; i < n_eq; ++i)
    {
      H_out[i] = hc::evaluate(F[i], point);
      for (size_t j = 0; j < var_idx_.size(); ++j)
        Hx_out[i][j] = hc::evaluate(J[i][var_idx_[j]], point);
      Complex ht(0);
      for (size_t k = 0; k < param_idx_.size(); ++k)
        ht += hc::evaluate(J[i][param_idx_[k]], point) * (p_target_[k] - p_start_[k]);
      Ht_out[i] = ht;
    }
  }

private:
  std::vector<Complex> full_point(const std::vector<Complex>& x, Scalar t) const
  {
    std::vector<Complex> point(F_.getVariables().size());
    for (size_t j = 0; j < var_idx_.size(); ++j)
      point[var_idx_[j]] = x[j];
    for (size_t k = 0; k < param_idx_.size(); ++k)
      point[param_idx_[k]] = (Scalar(1) - t) * p_start_[k] + t * p_target_[k];
    return point;
  }

  System<Scalar> F_;
  std::vector<int> var_idx_;
  std::vector<int> param_idx_;
  std::vector<Complex> p_start_;
  std::vector<Complex> p_target_;
};

} // namespace hc
