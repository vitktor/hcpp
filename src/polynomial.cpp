#include <hc/core/polynomial.hpp>

namespace hc
{
Polynomial<double> operator-(const Variable& var)
{
  return Polynomial<double>({-1.0}, {{1}}, {var});
}

Polynomial<double> operator+(const Variable& lv, const Variable& rv)
{
  if (lv == rv)
    return Polynomial<double>({2.0}, {{1}}, {lv});
  else if (lv < rv)
    return Polynomial<double>({1.0, 1.0}, {{1, 0}, {0, 1}}, {lv, rv});
  else
    return Polynomial<double>({1.0, 1.0}, {{1, 0}, {0, 1}}, {rv, lv});
}

Polynomial<double> operator-(const Variable& lv, const Variable& rv)
{
  if (lv == rv)
    return Polynomial<double>({}, {}, {lv});
  else if (lv < rv)
    return Polynomial<double>({1.0, -1.0}, {{1, 0}, {0, 1}}, {lv, rv});
  else
    return Polynomial<double>({-1.0, 1.0}, {{1, 0}, {0, 1}}, {rv, lv});
}

Polynomial<double> operator*(const Variable& lv, const Variable& rv)
{
  if (lv == rv)
    return Polynomial<double>({1.0}, {{2}}, {lv});
  else if (lv < rv)
    return Polynomial<double>({1.0}, {{1, 1}}, {lv, rv});
  else
    return Polynomial<double>({1.0}, {{1, 1}}, {rv, lv});
}

Polynomial<double> pow(const Variable& var, int exp)
{
  if (exp == 0)
    return Polynomial<double>(1.0);
  return Polynomial<double>({1.0}, {{exp}}, {var});
}
} // namespace hc
