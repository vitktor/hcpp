#pragma once
#include <complex>
#include <optional>

namespace hc
{

template <typename Scalar>
struct SolverOptions
{
  Scalar corrector_tol = Scalar(1e-12);
  int corrector_max_iters = 20;
  Scalar dt_init = Scalar(0.05);
  Scalar dt_min = Scalar(1e-9);
  Scalar dt_max = Scalar(0.1);
  int max_steps = 1000;
  Scalar max_norm = Scalar(1e6);
  bool verbose = false;
  // Worker threads for path tracking; 0 = one per hardware thread.
  int threads = 0;
  // Gamma for the straight-line homotopy; unset = a fresh random_gamma()
  // per solve. Fixing it makes repeated solves track identical paths.
  std::optional<std::complex<Scalar>> gamma;
};

} // namespace hc
