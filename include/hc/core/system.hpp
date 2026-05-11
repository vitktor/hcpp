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

private:
  std::vector<Polynomial<T>> polys;
  std::vector<Variable> vars;
};
} // namespace hc