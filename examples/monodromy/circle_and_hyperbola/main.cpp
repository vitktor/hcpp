#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <hc/hc.hpp>

using namespace hc;
using cd = std::complex<double>;

int main() {
  Variable x("x"), y("y"), a("a"), b("b");

  // Family F(x, y; a, b) = { x^2 + y^2 - a, x*y - b }: a circle intersected
  // with a hyperbola, 4 solutions for generic (a, b).
  System<double> F({pow(x, 2) + pow(y, 2) - a, x * y - b});
  std::vector<Variable> params{a, b};

  // Start pair: choose (x0, y0) freely, then (a0, b0) so that it solves F.
  cd x0(0.4, -0.2), y0(-0.9, 0.6);
  std::vector<cd> p0{x0 * x0 + y0 * y0, x0 * y0};

  MonodromySolver<double> solver;
  auto result = solver.solve(F, params, {x0, y0}, p0);

  std::ofstream out(std::string(HC_EXAMPLE_DIR) + "/results.txt");
  out << std::fixed << std::setprecision(6);
  out << "Monodromy on F = {x^2 + y^2 - a, x*y - b}\n\n";
  out << "Solutions at a = " << to_string(p0[0]) << ", b = " << to_string(p0[1])
      << " (" << result.loops << " loops, " << result.paths_tracked << " paths):\n";
  for (const auto& sol : result.solutions)
    out << "  x = " << to_string(sol[0]) << ", y = " << to_string(sol[1]) << "\n";

  // Move every solution to the parameters we actually care about.
  std::vector<cd> p1{cd(5.0), cd(2.0)};
  ParameterHomotopy<double> H(F, params, p0, p1);
  EulerPredictor<double> predictor;
  Tracker<double> tracker(H, predictor, 1e-12, 3, 0.05, 1e-9, 0.1, 1000);

  out << "\nSolutions at a = 5, b = 2:\n";
  for (const auto& sol : result.solutions) {
    auto r = tracker.track(sol, 0.0, 1.0);
    if (r.success)
      out << "  x = " << to_string(r.solution[0]) << ", y = " << to_string(r.solution[1]) << "\n";
  }

  std::cout << "Wrote results to " << HC_EXAMPLE_DIR << "/results.txt\n";
  return 0;
}
