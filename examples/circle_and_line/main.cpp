#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <hc/hc.hpp>

using namespace hc;
using cd = std::complex<double>;

// Manually predict+correct with a small fixed step, logging every step, so
// we can see exactly how the path behaves as it approaches wherever the
// full Tracker (adaptive step control) gives up.
static void trace_path(std::ofstream& out, const Homotopy<double>& H,
                       const Predictor<double>& predictor, std::vector<cd> point,
                       double t_end, double dt, int max_steps) {
  NewtonCorrector<double> corrector(1e-12, 20);
  double t = 0.0;
  for (int i = 0; i < max_steps && t < t_end; ++i) {
    std::vector<cd> predicted;
    predictor.predict(point, t, dt, H, predicted);
    double t_new = t + dt;
    auto result = corrector.correct(predicted, H, t_new);
    out << "  step " << i << ": t=" << t_new << "  x=" << predicted[0]
        << " y=" << predicted[1] << "  converged=" << result.converged
        << " iters=" << result.iterations << "\n";
    if (!result.converged) {
      out << "  -- Newton failed to converge here; stopping trace --\n";
      return;
    }
    point = predicted;
    t = t_new;
  }
}

int main() {
  Variable x("x"), y("y");

  // Start system G, degrees matching F's per-equation degrees (2 and 1) so
  // the Bezout number (2*1=2) matches F's actual solution count exactly:
  // G = { x^2 - 1, y - 1 },        roots (1,1), (-1,1)
  // F = { x^2 + y^2 - 1, x - y },  roots (+-1/sqrt(2), +-1/sqrt(2))
  auto g1 = pow(x, 2) - 1.0;
  auto g2 = y - 1.0;
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

  std::vector<std::vector<cd>> start_solutions{
      {cd(1.0, 0.0), cd(1.0, 0.0)},
      {cd(-1.0, 0.0), cd(1.0, 0.0)},
  };

  std::ofstream out(std::string(HC_EXAMPLE_DIR) + "/results.txt");
  out << std::fixed << std::setprecision(10);
  out << "Tracking G = {x^2 - 1, y - 1} to F = {x^2 + y^2 - 1, x - y}\n";
  out << "H(x,y,t) = (1-t)*G(x,y) + t*F(x,y), t: 0 -> 1\n\n";

  for (const auto& x0 : start_solutions) {
    auto result = tracker.track(x0, 0.0, 1.0);
    out << "start (x,y) = (" << x0[0] << ", " << x0[1] << ")\n";
    if (result.success) {
      out << "  -> solution (x,y) = (" << result.solution[0] << ", " << result.solution[1]
          << ")  (reached t=" << result.t << " in " << result.steps << " steps)\n\n";
    } else {
      out << "  -> FAILED (stalled at t=" << result.t << " after " << result.steps
          << " steps)\n\n";
    }
  }

  // Both paths above fail. Eliminating y from H2(x,y,t)=0 gives
  //   y = (1 - t - t*x) / (1 - 2*t)
  // which has a pole at t=1/2 unless x=1 exactly there. The trace below
  // confirms x drifts away from 1, so y genuinely diverges to infinity at
  // t=1/2: a real solution path leaving the finite plane, not a tolerance
  // artifact. This is the classical "path at infinity" failure of an
  // un-randomized straight-line homotopy; production HC software avoids it
  // with a random complex "gamma trick" (H = (1-t)*G + gamma*t*F), which
  // StraightLineHomotopy doesn't implement yet.
  out << "Manual trace from (1,1), dt=0.02:\n";
  trace_path(out, H, predictor, {cd(1.0, 0.0), cd(1.0, 0.0)}, 0.5, 0.02, 60);

  // Gamma trick: H = gamma*(1-t)*G + t*F for a fixed non-real unit gamma
  // steers the paths off the real axis, avoiding the singularities above
  // (see StraightLineHomotopy's gamma parameter).
  cd gamma(0.6, 0.8);
  StraightLineHomotopy<double> H_gamma(start, target, {.gamma = gamma});
  Tracker<double> tracker_gamma(H_gamma, predictor, 1e-12, 20, 0.05, 1e-9, 0.1, 1000);
  out << "\nWith gamma = " << gamma << ":\n";
  for (const auto& x0 : start_solutions) {
    auto result = tracker_gamma.track(x0, 0.0, 1.0);
    out << "start (x,y) = (" << x0[0] << ", " << x0[1] << ")\n";
    if (result.success) {
      out << "  -> solution (x,y) = (" << result.solution[0] << ", " << result.solution[1]
          << ")  (reached t=" << result.t << " in " << result.steps << " steps)\n\n";
    } else {
      out << "  -> FAILED (stalled at t=" << result.t << " after " << result.steps
          << " steps)\n\n";
    }
  }

  std::cout << "Wrote results to " << HC_EXAMPLE_DIR << "/results.txt\n";
  return 0;
}
