#pragma once

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
};

} // namespace hc
