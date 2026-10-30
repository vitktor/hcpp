# hcpp
A C++ library for solving systems of polynomial equations using numerical homotopy continuation.

## Solver

Solve `F = {x^2 + y^2 - 1, x^2 - y}` (a circle intersected with a parabola).
The Bezout number is 4, but only 2 of the 4 solutions are real. `Solver`
finds all of them:

```cpp
#include <hc/hc.hpp>
#include <iostream>

using namespace hc;

int main() {
  Variable x("x"), y("y");
  System<double> F({pow(x, 2) + pow(y, 2) - 1.0, pow(x, 2) - y}, {x, y});

  Solver<double> solver;
  auto result = solver.solve(F);

  for (const auto& sol : result.solutions())
    std::cout << to_string(sol[0]) << ", " << to_string(sol[1]) << "\n";
  // 0.786151+0.000000j, 0.618034+0.000000j
  // 0.000000+1.272020j, -1.618034+0.000000j
  // -0.786151-0.000000j, 0.618034+0.000000j
  // 0.000000-1.272020j, -1.618034+0.000000j
}
```

See `examples/solving` for more, including filtering to just the real
solutions via `result.real_solutions()`.

## Tracking

`Solver` builds its start system and tracks every path automatically. You
can also drive a single path yourself for more control. A 2-equation,
2-variable system: track `G = {x^2 - 1, y^2 - 1}`
(roots `(+-1,+-1)`) to `F = {x^2 + y^2 - 1, x - y}` (roots
`(+-1/sqrt(2), +-1/sqrt(2))`). `G` has 4 roots but `F` only has 2, so two of
the four paths necessarily diverge to infinity. We track only the two that converge, `(1,1)` and `(-1,-1)`:

```cpp
#include <hc/hc.hpp>
#include <iostream>

using namespace hc;
using cd = std::complex<double>;

int main() {
  Variable x("x"), y("y");

  auto g1 = pow(x, 2) - 1.0;
  auto g2 = pow(y, 2) - 1.0;
  auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
  auto f2 = x - y;

  System<double> start({g1, g2}, {x, y});
  System<double> target({f1, f2}, {x, y});
  StraightLineHomotopy<double> H(start, target);
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor,
                          /*corrector_tol=*/1e-12, /*corrector_max_iters=*/20,
                          /*dt_init=*/0.05, /*dt_min=*/1e-9, /*dt_max=*/0.1,
                          /*max_steps=*/1000);

  for (auto x0 : std::vector<std::vector<cd>>{{cd(1, 0), cd(1, 0)}, {cd(-1, 0), cd(-1, 0)}}) {
    auto result = tracker.track(x0, 0.0, 1.0);
    std::cout << to_string(result.solution[0]) << ", " << to_string(result.solution[1]) << "\n";
  }
  // 0.707107+0.000000j, 0.707107+0.000000j
  // -0.707107-0.000000j, -0.707107-0.000000j
}
```

See `examples/tracking` for more, including the "gamma trick" that
`StraightLineHomotopy` uses by default to keep paths away from singularities.

## Monodromy

For a parametric family `F(x; p)`, `MonodromySolver` finds all solutions at a
generic parameter point from just one known solution: it tracks solutions
around random loops in parameter space, and since the loops permute the
solutions, every endpoint is again a solution. `ParameterHomotopy` then moves
them to the parameters you actually care about. Parameters are ordinary
`Variable`s of the `System`; you only say which ones they are:

```cpp
#include <hc/hc.hpp>
#include <iostream>

using namespace hc;
using cd = std::complex<double>;

int main() {
  Variable x("x"), y("y"), a("a"), b("b");
  System<double> F({pow(x, 2) + pow(y, 2) - a, x * y - b});
  std::vector<Variable> params{a, b};

  // Start pair: choose (x0, y0) freely, then (a0, b0) so that it solves F.
  cd x0(0.4, -0.2), y0(-0.9, 0.6);
  std::vector<cd> p0{x0 * x0 + y0 * y0, x0 * y0};

  MonodromySolver<double> solver;
  auto result = solver.solve(F, params, {x0, y0}, p0);  // 4 solutions at p0

  ParameterHomotopy<double> H(F, params, p0, {cd(5.0), cd(2.0)});
  EulerPredictor<double> predictor;
  Tracker<double> tracker(H, predictor, 1e-12, 3, 0.05, 1e-9, 0.1, 1000);
  for (const auto& sol : result.solutions) {
    auto r = tracker.track(sol, 0.0, 1.0);
    std::cout << to_string(r.solution[0]) << ", " << to_string(r.solution[1]) << "\n";
  }
  // -1.000000+0.000000j, -2.000000-0.000000j
  // 2.000000+0.000000j, 1.000000-0.000000j
  // -2.000000-0.000000j, -1.000000+0.000000j
  // 1.000000-0.000000j, 2.000000+0.000000j
}
```

Without `target_solutions_count`, the solver stops once
`max_loops_no_progress` new loops in a row found nothing new, so the result
can be incomplete; set the count when you know it. See `examples/monodromy`.
