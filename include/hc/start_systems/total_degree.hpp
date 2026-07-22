#pragma once
#include <cmath>
#include <complex>
#include <numbers>
#include <utility>
#include <vector>

#include <hc/core/polynomial.hpp>
#include <hc/core/system.hpp>
#include <hc/start_systems/start_system.hpp>

namespace hc
{

// G_i(x) = x_i^{d_i} - 1, where d_i = deg(target_i). Start solutions are
// the Cartesian product of each equation's d_i-th roots of unity: the full
// Bezout number (product of degrees) of start points. Track these through
// a StraightLineHomotopy(total_degree_system.system(), target); its
// default random gamma provides the genericity needed to avoid path
// crossings and paths at infinity for any excess (Bezout minus actual
// root count) paths.
template <typename Scalar>
class TotalDegreeSystem : public StartSystem<Scalar>
{
public:
  explicit TotalDegreeSystem(const System<Scalar>& target)
      : system_(build_system(target))
  {
    build_start_solutions(target);
  }

  const System<Scalar>& system() const override { return system_; }

  const std::vector<std::vector<std::complex<Scalar>>>& start_solutions() const override
  {
    return start_solutions_;
  }

private:
  static System<Scalar> build_system(const System<Scalar>& target)
  {
    const auto& vars = target.getVariables();
    const auto& polys = target.getPolynomials();
    std::vector<Polynomial<Scalar>> g_polys;
    g_polys.reserve(polys.size());
    for (size_t i = 0; i < polys.size(); ++i)
      g_polys.push_back(pow(vars[i], polys[i].degree()) - Scalar(1));
    return System<Scalar>(std::move(g_polys), vars);
  }

  void build_start_solutions(const System<Scalar>& target)
  {
    start_solutions_ = {{}};
    for (const auto& p : target.getPolynomials())
    {
      int d = p.degree();
      std::vector<std::complex<Scalar>> roots;
      roots.reserve(d);
      for (int k = 0; k < d; ++k)
      {
        Scalar theta = Scalar(2) * std::numbers::pi_v<Scalar> * Scalar(k) / Scalar(d);
        roots.emplace_back(std::cos(theta), std::sin(theta));
      }

      std::vector<std::vector<std::complex<Scalar>>> next;
      next.reserve(start_solutions_.size() * roots.size());
      for (const auto& partial : start_solutions_)
        for (const auto& r : roots)
        {
          auto sol = partial;
          sol.push_back(r);
          next.push_back(std::move(sol));
        }
      start_solutions_ = std::move(next);
    }
  }

  System<Scalar> system_;
  std::vector<std::vector<std::complex<Scalar>>> start_solutions_;
};

} // namespace hc
