#pragma once
#include <complex>
#include <vector>

#include <hc/core/system.hpp>

namespace hc
{

template <typename Scalar>
class StartSystem
{
public:
  virtual ~StartSystem() = default;
  virtual const System<Scalar>& system() const = 0;
  virtual const std::vector<std::vector<std::complex<Scalar>>>& start_solutions() const = 0;
};

} // namespace hc
