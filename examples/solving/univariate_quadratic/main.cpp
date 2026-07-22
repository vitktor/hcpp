#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <hc/hc.hpp>

using namespace hc;

int main() {
  Variable x("x");

  // F(x) = x^2 - 9, roots +-3. Unlike examples/tracking, we don't build a
  // start system or specify start points by hand -- Solver builds a
  // TotalDegreeSystem and tracks every one of its start solutions itself.
  System<double> target({pow(x, 2) - 9.0}, {x});

  Solver<double> solver;
  auto result = solver.solve(target);

  std::ofstream out(std::string(HC_EXAMPLE_DIR) + "/results.txt");
  out << std::fixed << std::setprecision(6);
  out << "Solving F(x) = x^2 - 9 = 0\n\n";

  for (const auto& sol : result.real_solutions())
    out << "x = " << to_string(sol[0]) << "\n";

  std::cout << "Wrote results to " << HC_EXAMPLE_DIR << "/results.txt\n";
  return 0;
}
