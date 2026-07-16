#include <hc/core/polynomial.hpp>

namespace hc
{
Variable::operator Polynomial<double>() const
{
  return Polynomial<double>({1.0}, {{1}}, {*this});
}

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

Polynomial<double> operator+(const Polynomial<double>& poly, double scalar)
{
  return poly + Polynomial<double>(scalar);
}

Polynomial<double> operator+(double scalar, const Polynomial<double>& poly)
{
  return poly + Polynomial<double>(scalar);
}

Polynomial<double> operator-(const Polynomial<double>& poly, double scalar)
{
  return poly - Polynomial<double>(scalar);
}

Polynomial<double> operator-(double scalar, const Polynomial<double>& poly)
{
  return Polynomial<double>(scalar) - poly;
}

Polynomial<double> operator*(const Polynomial<double>& poly, double scalar)
{
  return poly * Polynomial<double>(scalar);
}

Polynomial<double> operator*(double scalar, const Polynomial<double>& poly)
{
  return poly * Polynomial<double>(scalar);
}
} // namespace hc
