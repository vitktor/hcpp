#pragma once
#include <algorithm>
#include <vector>

#include <hc/core/variable.hpp>

namespace hc
{
template <typename T> class Polynomial
{
public:
  Polynomial(std::vector<T> coefficients,
             std::vector<std::vector<int>> exponents,
             std::vector<Variable> variables, bool sorted_terms = false)
      : coeffs(std::move(coefficients)),
        exps(std::move(exponents)),
        vars(std::move(variables))
  {
    if (!sorted_terms)
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
              [&](size_t a, size_t b) { return exps[a] > exps[b]; });
    std::vector<std::vector<int>> sorted_exps(n);
    std::vector<T> sorted_coeffs(n);
    for (size_t i = 0; i < n; ++i)
    {
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
    for (const auto& e : exps)
    {
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

    std::vector<std::vector<int>> result_exps;
    std::vector<T> result_coeffs;
    result_exps.reserve(exps.size() + poly.exps.size());
    result_coeffs.reserve(coeffs.size() + poly.coeffs.size());

    size_t i = 0, j = 0;
    size_t n1 = exps.size();
    size_t n2 = poly.exps.size();

    while (i < n1 && j < n2)
    {
      auto ki = map_exp(exps[i], mvars.map1, n_vars);
      auto kj = map_exp(poly.exps[j], mvars.map2, n_vars);
      if (kj < ki)
      {
        result_exps.push_back(std::move(ki));
        result_coeffs.push_back(coeffs[i++]);
      }
      else if (ki < kj)
      {
        result_exps.push_back(std::move(kj));
        result_coeffs.push_back(poly.coeffs[j++]);
      }
      else
      {
        T sum = coeffs[i++] + poly.coeffs[j++];
        if (sum != T(0))
        {
          result_exps.push_back(std::move(ki));
          result_coeffs.push_back(sum);
        }
      }
    }
    while (i < n1)
    {
      result_exps.push_back(map_exp(exps[i], mvars.map1, n_vars));
      result_coeffs.push_back(coeffs[i++]);
    }
    while (j < n2)
    {
      result_exps.push_back(map_exp(poly.exps[j], mvars.map2, n_vars));
      result_coeffs.push_back(poly.coeffs[j++]);
    }

    constexpr bool sorted_terms = true;
    return Polynomial<T>(std::move(result_coeffs), std::move(result_exps),
                         std::move(mvars.vars), sorted_terms);
  }

  Polynomial<T> operator-(const Polynomial<T>& poly) const
  {
    return *this + (-poly);
  }

  Polynomial<T> operator*(const Polynomial<T>& poly) const
  {
    MergedVars mvars = merge_sorted_vars(vars, poly.vars);
    size_t n_vars = mvars.vars.size();

    size_t n1 = exps.size(), n2 = poly.exps.size();
    std::vector<std::vector<int>> prod_exps;
    std::vector<T> prod_coeffs;
    prod_exps.reserve(n1 * n2);
    prod_coeffs.reserve(n1 * n2);

    for (size_t i = 0; i < n1; ++i)
    {
      auto ki = map_exp(exps[i], mvars.map1, n_vars);
      for (size_t j = 0; j < n2; ++j)
      {
        auto kj = map_exp(poly.exps[j], mvars.map2, n_vars);
        std::vector<int> e(n_vars);
        for (size_t k = 0; k < n_vars; ++k)
          e[k] = ki[k] + kj[k];
        prod_exps.push_back(std::move(e));
        prod_coeffs.push_back(coeffs[i] * poly.coeffs[j]);
      }
    }

    size_t n = prod_exps.size();
    std::vector<size_t> idx(n);
    for (size_t k = 0; k < n; ++k)
      idx[k] = k;
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return prod_exps[a] > prod_exps[b]; });

    std::vector<std::vector<int>> result_exps;
    std::vector<T> result_coeffs;
    result_exps.reserve(n);
    result_coeffs.reserve(n);

    for (size_t k = 0; k < n; ++k)
    {
      auto& e = prod_exps[idx[k]];
      T c = prod_coeffs[idx[k]];
      if (!result_exps.empty() && result_exps.back() == e)
      {
        result_coeffs.back() += c;
      }
      else
      {
        if (!result_exps.empty() && result_coeffs.back() == T(0))
        {
          result_exps.pop_back();
          result_coeffs.pop_back();
        }
        result_exps.push_back(std::move(e));
        result_coeffs.push_back(c);
      }
    }
    if (!result_exps.empty() && result_coeffs.back() == T(0))
    {
      result_exps.pop_back();
      result_coeffs.pop_back();
    }

    return Polynomial<T>(std::move(result_coeffs), std::move(result_exps),
                         std::move(mvars.vars), true);
  }

private:
  std::vector<T> coeffs;
  std::vector<std::vector<int>> exps;
  std::vector<Variable> vars;

  static std::vector<int> map_exp(const std::vector<int>& e,
                                  const std::vector<int>& map, size_t n_vars)
  {
    std::vector<int> key(n_vars, 0);
    for (size_t k = 0; k < map.size(); ++k)
      key[map[k]] = e[k];
    return key;
  }
};

Polynomial<double> operator-(const Variable& var);
Polynomial<double> operator+(const Variable& lv, const Variable& rv);
Polynomial<double> operator-(const Variable& lv, const Variable& rv);
Polynomial<double> operator*(const Variable& lv, const Variable& rv);
Polynomial<double> pow(const Variable& var, int exp);

} // namespace hc
