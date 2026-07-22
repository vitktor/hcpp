#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <hc/hc.hpp>

using namespace hc;

int main() {
  Variable x("x"), y("y");

  // F = { x^2 + y^2 - 1, x^2 - y }: a circle intersected with a parabola.
  // Bezout number 2*2=4, but only 2 of the 4 solutions are real -- the
  // other 2 have a genuinely complex x. Solver finds all 4; real_solutions()
  // filters to the 2 that are actually real.
  System<double> target({pow(x, 2) + pow(y, 2) - 1.0, pow(x, 2) - y}, {x, y});

  Solver<double> solver;
  auto result = solver.solve(target);

  std::ofstream out(std::string(HC_EXAMPLE_DIR) + "/results.txt");
  out << std::fixed << std::setprecision(6);
  out << "Solving F = {x^2 + y^2 - 1, x^2 - y} = 0\n\n";

  out << "All solutions:\n";
  for (const auto& sol : result.solutions())
    out << "  x = " << to_string(sol[0]) << ", y = " << to_string(sol[1]) << "\n";

  out << "\nReal solutions:\n";
  for (const auto& sol : result.real_solutions())
    out << "  x = " << to_string(sol[0]) << ", y = " << to_string(sol[1]) << "\n";

  std::cout << "Wrote results to " << HC_EXAMPLE_DIR << "/results.txt\n";
  return 0;
}
