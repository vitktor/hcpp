#pragma once
#include <algorithm>
#include <complex>
#include <vector>

#include <hc/tracking/tracker.hpp>

namespace hc
{

// Ties every start solution's TrackResult together. Doesn't deduplicate
// solutions reached by more than one path (a real possibility whenever the
// start system's Bezout number exceeds the target's true root count), and
// doesn't distinguish regular from singular solutions -- that needs an
// endgame (condition number / winding number), not implemented yet.
template <typename Scalar>
struct Result
{
  std::vector<TrackResult<Scalar>> paths;

  std::vector<std::vector<std::complex<Scalar>>> solutions() const
  {
    std::vector<std::vector<std::complex<Scalar>>> sols;
    for (const auto& p : paths)
      if (p.success)
        sols.push_back(p.solution);
    return sols;
  }

  std::vector<std::vector<std::complex<Scalar>>> real_solutions(Scalar tol = Scalar(1e-8)) const
  {
    std::vector<std::vector<std::complex<Scalar>>> sols;
    for (const auto& p : paths)
    {
      if (!p.success)
        continue;
      bool real = std::all_of(p.solution.begin(), p.solution.end(),
                              [tol](const std::complex<Scalar>& xi) {
                                return std::abs(xi.imag()) < tol;
                              });
      if (real)
        sols.push_back(p.solution);
    }
    return sols;
  }
};

} // namespace hc
