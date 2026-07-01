#pragma once
#include <cmath>
#include <complex>

#include <hc/core/polynomial.hpp>

namespace hc
{

template <typename T>
Polynomial<T> differentiate(const Polynomial<T>& poly, const Variable& var)
{
  const auto& vars = poly.getVariables();
  auto it = std::find(vars.begin(), vars.end(), var);
  if (it == vars.end())
    return Polynomial<T>({}, {}, vars, true);

  int var_idx = it - vars.begin();
  const auto& exps = poly.getExponents();
  const auto& coeffs = poly.getCoefficients();

  std::vector<std::vector<int>> result_exps;
  std::vector<T> result_coeffs;

  for (size_t i = 0; i < exps.size(); ++i) {
    int e = exps[i][var_idx];
    if (e == 0)
      continue;
    auto new_exp = exps[i];
    new_exp[var_idx] -= 1;
    result_exps.push_back(std::move(new_exp));
    result_coeffs.push_back(coeffs[i] * T(e));
  }

  if (result_exps.empty())
    return Polynomial<T>({}, {}, vars, true);

  return Polynomial<T>(std::move(result_coeffs), std::move(result_exps), vars, true);
}

template <typename T, typename S>
std::common_type_t<T, S> evaluate(const Polynomial<T>& poly, const std::vector<S>& point)
{
  using R = std::common_type_t<T, S>;
  const auto& exps = poly.getExponents();
  const auto& coeffs = poly.getCoefficients();

  R result = R(0);
  for (size_t i = 0; i < coeffs.size(); ++i) {
    R term = static_cast<R>(coeffs[i]);
    for (size_t j = 0; j < point.size(); ++j)
      term *= std::pow(point[j], exps[i][j]);
    result += term;
  }
  return result;
}

} // namespace hc
