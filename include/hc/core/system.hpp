#pragma once
#include <vector>

#include <hc/core/calculus.hpp>
#include <hc/core/polynomial.hpp>

namespace hc
{
template <typename T> class System
{
public:
  System(std::vector<Polynomial<T>> polys, std::vector<Variable> vars)
      : polys(std::move(polys)), vars(std::move(vars))
  {
    this->polys = update_vars(std::move(this->polys), this->vars);
    jac = compute_jacobian();
  }

  System(std::vector<Polynomial<T>> polys) : polys(std::move(polys))
  {
    for (const auto& p : this->polys)
      vars.insert(vars.end(), p.getVariables().begin(), p.getVariables().end());
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
    this->polys = update_vars(std::move(this->polys), vars);
    jac = compute_jacobian();
  }

  const std::vector<Polynomial<T>>& getPolynomials() const { return polys; }
  const std::vector<Variable>& getVariables() const { return vars; }
  const std::vector<std::vector<Polynomial<T>>>& getJacobian() const { return jac; }

private:
  std::vector<Polynomial<T>> polys;
  std::vector<Variable> vars;
  std::vector<std::vector<Polynomial<T>>> jac;

  std::vector<std::vector<Polynomial<T>>> compute_jacobian() const
  {
    std::vector<std::vector<Polynomial<T>>> J(
        polys.size(), std::vector<Polynomial<T>>(vars.size(), Polynomial<T>(T(0))));
    for (size_t i = 0; i < polys.size(); ++i)
      for (size_t j = 0; j < vars.size(); ++j)
        J[i][j] = differentiate(polys[i], vars[j]);
    return J;
  }
};

} // namespace hc