#pragma once
#include <vector>

#include <hc/core/polynomial.hpp>

namespace hc
{
template <typename T> class System
{
public:
  System(std::vector<Polynomial<T>> polys, std::vector<Variable> vars)
      : polys(std::move(polys)), vars(std::move(vars))
  {
  }

  System(std::vector<Polynomial<T>> polys) : polys(std::move(polys))
  {
    for (const auto& p : this->polys)
      vars.insert(vars.end(), p.getVariables().begin(), p.getVariables().end());
    std::sort(vars.begin(), vars.end());
    vars.erase(std::unique(vars.begin(), vars.end()), vars.end());
  }

  const std::vector<Polynomial<T>>& getPolynomials() const { return polys; }
  const std::vector<Variable>& getVariables() const { return vars; }

private:
  std::vector<Polynomial<T>> polys;
  std::vector<Variable> vars;
};
} // namespace hc