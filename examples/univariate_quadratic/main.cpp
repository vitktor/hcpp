#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <hc/hc.hpp>

using namespace hc;
using cd = std::complex<double>;

int main() {
  Variable x("x");

  // G(x) = x^2 - 4, roots {2, -2}
  // F(x) = x^2 - 9, roots {3, -3}
  auto g = pow(x, 2) - 4.0;
  auto f = pow(x, 2) - 9.0;

  System<double> start({g}, {x});
  System<double> target({f}, {x});
  StraightLineHomotopy<double> H(start, target);
  EulerPredictor<double> predictor;

  Tracker<double> tracker(H, predictor,
                          /*corrector_tol=*/1e-12, /*corrector_max_iters=*/20,
                          /*dt_init=*/0.05, /*dt_min=*/1e-9, /*dt_max=*/0.1,
                          /*max_steps=*/1000);

  std::vector<cd> start_solutions{cd(2.0, 0.0), cd(-2.0, 0.0)};

  std::ofstream out(std::string(HC_EXAMPLE_DIR) + "/results.txt");
  out << std::fixed << std::setprecision(10);
  out << "Tracking G(x) = x^2 - 4 (roots +-2) to F(x) = x^2 - 9 (roots +-3)\n";
  out << "H(x,t) = (1-t)*G(x) + t*F(x), t: 0 -> 1\n\n";

  for (const auto& x0 : start_solutions) {
    auto result = tracker.track({x0}, 0.0, 1.0);
    out << "start x = " << x0 << "\n";
    if (result.success) {
      out << "  -> solution x = " << result.solution[0]
          << "  (reached t=" << result.t << " in " << result.steps << " steps)\n\n";
    } else {
      out << "  -> FAILED (stalled at t=" << result.t << " after " << result.steps
          << " steps)\n\n";
    }
  }

  std::cout << "Wrote results to " << HC_EXAMPLE_DIR << "/results.txt\n";
  return 0;
}
