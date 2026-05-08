#pragma once
#include <algorithm>
#include <vector>

#include <hc/core/variable.hpp>

namespace hc
{
template <typename T> class Polynomial
{
public:
  Polynomial(std::vector<T> coefficients, std::vector<std::vector<int>> exponents,
             std::vector<Variable> variables, bool presorted = false)
      : coeffs(std::move(coefficients)), exps(std::move(exponents)),
        vars(std::move(variables))
  {
    if (!presorted)
      sort_terms();
  }

  explicit Polynomial(T scalar) : coeffs({scalar}), exps({{}}), vars({}) {}

  void sort_terms()
  {
    size_t n = exps.size();
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i)
      idx[i] = i;
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return exps[a] < exps[b]; });
    std::vector<std::vector<int>> sorted_exps(n);
    std::vector<T> sorted_coeffs(n);
    for (size_t i = 0; i < n; ++i) {
      sorted_exps[i] = std::move(exps[idx[i]]);
      sorted_coeffs[i] = std::move(coeffs[idx[i]]);
    }
    exps = std::move(sorted_exps);
    coeffs = std::move(sorted_coeffs);
  }

  const std::vector<T>& getCoefficients() const { return coeffs; }

  void setCoefficients(std::vector<T> coefficients)
  {
    coeffs = std::move(coefficients);
  }

  const std::vector<std::vector<int>>& getExponents() const { return exps; }

  void setExponents(std::vector<std::vector<int>> exponents)
  {
    exps = std::move(exponents);
  }

  const std::vector<Variable>& getVariables() const { return vars; }

  void setVariables(std::vector<Variable> variables)
  {
    vars = std::move(variables);
  }

  int degree() const
  {
    if (exps.empty())
      return -1;
    int max_deg = 0;
    for (const auto& e : exps) {
      int deg = 0;
      for (int x : e)
        deg += x;
      max_deg = std::max(max_deg, deg);
    }
    return max_deg;
  }

  int degree(const Variable& var) const
  {
    auto it = std::find(vars.begin(), vars.end(), var);
    if (it == vars.end())
      return 0;
    int id = it - vars.begin();
    int max_exp = 0;
    for (const auto& e : exps)
      max_exp = std::max(max_exp, e[id]);
    return max_exp;
  }

  Polynomial<T> operator-() const
  {
    std::vector<T> neg_coeffs(coeffs.size());
    std::transform(coeffs.begin(), coeffs.end(), neg_coeffs.begin(),
                   [](const T& c) { return -c; });
    return Polynomial<T>(std::move(neg_coeffs), exps, vars);
  }

  Polynomial<T> operator+(const Polynomial<T>& poly) const
  {
    MergedVars mvars = merge_sorted_vars(getVariables(), poly.getVariables());
    size_t n_vars = mvars.vars.size();

    auto remap = [&](const std::vector<int>& e, const std::vector<int>& map) {
      std::vector<int> key(n_vars, 0);
      for (size_t i = 0; i < map.size(); ++i)
        key[map[i]] = e[i];
      return key;
    };

    std::vector<std::vector<int>> result_exps;
    std::vector<T> result_coeffs;
    result_exps.reserve(exps.size() + poly.exps.size());
    result_coeffs.reserve(coeffs.size() + poly.coeffs.size());

    size_t i = 0, j = 0;
    size_t n1 = exps.size();
    size_t n2 = poly.exps.size();

    while (i < n1 && j < n2) {
      auto ki = remap(exps[i], mvars.map1);
      auto kj = remap(poly.exps[j], mvars.map2);
      if (ki < kj) {
        result_exps.push_back(std::move(ki));
        result_coeffs.push_back(coeffs[i++]);
      } else if (kj < ki) {
        result_exps.push_back(std::move(kj));
        result_coeffs.push_back(poly.coeffs[j++]);
      } else {
        T sum = coeffs[i++] + poly.coeffs[j++];
        if (sum != T(0)) {
          result_exps.push_back(std::move(ki));
          result_coeffs.push_back(sum);
        }
      }
    }
    while (i < n1) {
      result_exps.push_back(remap(exps[i], mvars.map1));
      result_coeffs.push_back(coeffs[i++]);
    }
    while (j < n2) {
      result_exps.push_back(remap(poly.exps[j], mvars.map2));
      result_coeffs.push_back(poly.coeffs[j++]);
    }

    if (result_exps.empty())
      return Polynomial<T>(0);

    return Polynomial<T>(std::move(result_coeffs), std::move(result_exps),
                         std::move(mvars.vars), true);
  }

  Polynomial<T> operator-(const Polynomial<T>& poly) const
  {
  }

  Polynomial<T> operator*(const Polynomial<T>& poly) const
  {
  }

private:
  std::vector<T> coeffs;
  std::vector<std::vector<int>> exps;
  std::vector<Variable> vars;
};

Polynomial<double> operator-(const Variable& var);
Polynomial<double> operator+(const Variable& lv, const Variable& rv);
Polynomial<double> operator-(const Variable& lv, const Variable& rv);
Polynomial<double> operator*(const Variable& lv, const Variable& rv);
Polynomial<double> pow(const Variable& var, int exp);

} // namespace hc
