#include <algorithm>
#include <atomic>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>
#include <hc/parallel.hpp>
#include <hc/solver/monodromy.hpp>
#include <hc/solver/solver.hpp>

using namespace hc;
using cd = std::complex<double>;

TEST(ParallelFor, VisitsEveryIndexOnce) {
  for (int threads : {1, 2, 4, 8, 0}) {
    for (size_t n : {0u, 1u, 3u, 1000u}) {
      std::vector<std::atomic<int>> hits(n);
      std::atomic<int> max_tid{0};
      parallel_for(n, threads, [&](int tid, size_t i) {
        hits[i].fetch_add(1);
        int m = max_tid.load();
        while (tid > m && !max_tid.compare_exchange_weak(m, tid)) {}
      });
      for (size_t i = 0; i < n; ++i)
        EXPECT_EQ(hits[i].load(), 1) << "threads=" << threads << " n=" << n << " i=" << i;
      EXPECT_LT(max_tid.load(), std::max(1, resolve_threads(threads)));
    }
  }
}

TEST(ParallelFor, RethrowsWorkerException) {
  EXPECT_THROW(parallel_for(100, 4,
                            [](int, size_t i) {
                              if (i == 37)
                                throw std::runtime_error("boom");
                            }),
               std::runtime_error);
}

namespace {

// 5 coupled quadratics: Bezout number 2^5 = 32 paths.
System<double> coupled_quadratics() {
  std::vector<Variable> x;
  for (int i = 0; i < 5; ++i)
    x.emplace_back("x" + std::to_string(i));
  std::vector<Polynomial<double>> eqs;
  for (int i = 0; i < 5; ++i)
    eqs.push_back(pow(x[i], 2) + 0.5 * (x[i] * x[(i + 1) % 5]) - 0.3 * Polynomial<double>(x[(i + 2) % 5]) - 1.0);
  return System<double>(eqs, x);
}

} // namespace

TEST(ParallelSolver, ThreadCountDoesNotChangeResult) {
  // With gamma fixed, every thread count must track exactly the same paths
  // and return them in start-solution order.
  auto F = coupled_quadratics();
  SolverOptions<double> opts;
  opts.gamma = cd(0.6, 0.8);
  opts.threads = 1;
  auto serial = Solver<double>(opts).solve(F);
  ASSERT_EQ(serial.paths.size(), 32u);
  EXPECT_GT(serial.solutions().size(), 0u);

  for (int threads : {2, 8, 0}) {
    opts.threads = threads;
    auto par = Solver<double>(opts).solve(F);
    ASSERT_EQ(par.paths.size(), serial.paths.size());
    for (size_t i = 0; i < serial.paths.size(); ++i) {
      EXPECT_EQ(par.paths[i].success, serial.paths[i].success) << "path " << i;
      EXPECT_EQ(par.paths[i].status, serial.paths[i].status) << "path " << i;
      EXPECT_EQ(par.paths[i].steps, serial.paths[i].steps) << "path " << i;
      EXPECT_EQ(par.paths[i].solution, serial.paths[i].solution) << "path " << i; // bitwise
    }
  }
}

TEST(ParallelSolver, VerboseProgressLinesAreIntact) {
  // 32 paths on 8 threads all report to one stream. Every line must be
  // whole (no interleaving), the finished counter must run 1..32 in order,
  // and every path must be reported exactly once.
  auto F = coupled_quadratics();
  std::ostringstream log;
  SolverOptions<double> opts;
  opts.threads = 8;
  opts.verbose = true;
  opts.log = &log;
  auto result = Solver<double>(opts).solve(F);

  std::regex line_re(R"(\[(\d+)/32\] path (\d+): (success|diverged|failed), (\d+) steps)");
  std::istringstream lines(log.str());
  std::string line;
  size_t count = 0;
  std::set<size_t> paths_seen;
  while (std::getline(lines, line)) {
    std::smatch m;
    ASSERT_TRUE(std::regex_match(line, m, line_re)) << "garbled line: " << line;
    ++count;
    EXPECT_EQ(std::stoul(m[1]), count);
    size_t path = std::stoul(m[2]);
    EXPECT_TRUE(paths_seen.insert(path).second) << "path " << path << " reported twice";
    EXPECT_EQ(m[3], to_string(result.paths[path].status));
    EXPECT_EQ(std::stoi(m[4]), result.paths[path].steps);
  }
  EXPECT_EQ(count, 32u);
}

TEST(ParallelSolver, QuietByDefault) {
  auto F = coupled_quadratics();
  std::ostringstream log;
  SolverOptions<double> opts;
  opts.log = &log;
  Solver<double>(opts).solve(F);
  EXPECT_TRUE(log.str().empty());
}

TEST(ParallelMonodromy, ThreadCountDoesNotChangeResult) {
  // x^6 + a*x + b: 6 solutions for generic (a, b).
  Variable x("x"), a("a"), b("b");
  System<double> F({pow(x, 6) + a * x + b});
  std::vector<Variable> params{a, b};
  cd x0(0.3, 0.7), a0(-0.5, 1.1);
  cd b0 = -(std::pow(x0, 6) + a0 * x0);

  MonodromyOptions<double> opts;
  opts.threads = 1;
  auto serial = MonodromySolver<double>(opts).solve(F, params, {x0}, {a0, b0});
  EXPECT_EQ(serial.solutions.size(), 6u);

  for (int threads : {2, 8}) {
    opts.threads = threads;
    auto par = MonodromySolver<double>(opts).solve(F, params, {x0}, {a0, b0});
    EXPECT_EQ(par.loops, serial.loops);
    EXPECT_EQ(par.paths_tracked, serial.paths_tracked);
    ASSERT_EQ(par.solutions.size(), serial.solutions.size());
    for (size_t i = 0; i < serial.solutions.size(); ++i)
      EXPECT_EQ(par.solutions[i], serial.solutions[i]); // bitwise identical
  }
}
