#pragma once
#include <algorithm>
#include <complex>
#include <cstdint>
#include <random>
#include <vector>

#include <hc/core/system.hpp>
#include <hc/homotopy/parameter.hpp>
#include <hc/tracking/runge_kutta.hpp>
#include <hc/tracking/tracker.hpp>

namespace hc
{

template <typename Scalar>
struct MonodromyOptions
{
  Scalar corrector_tol = Scalar(1e-12);
  // Few iterations on purpose: a step is only accepted if Newton converges
  // quickly, which keeps the tracker from jumping onto a neighboring path.
  // With 20 iterations, loops of x^6 + a*x + b routinely sent two roots to
  // the same endpoint, so loops no longer acted as permutations.
  int corrector_max_iters = 3;
  Scalar dt_init = Scalar(0.05);
  Scalar dt_min = Scalar(1e-9);
  Scalar dt_max = Scalar(0.1);
  int max_steps = 1000;
  Scalar max_norm = Scalar(1e6);

  // Stop as soon as this many solutions are known (0 = unknown, rely on
  // max_loops_no_progress alone).
  int target_solutions_count = 0;
  // Stop after this many consecutive new loops that produced no new solution.
  // Higher than HC.jl's 5: a random loop often fixes a given solution, and
  // with 5, x^6 + a*x + b and x^10 + a*x + b stopped early for ~1 in 4 seeds.
  int max_loops_no_progress = 10;
  // Hard cap on the total number of loops.
  int max_loops = 100;
  // Two endpoints closer than this (max-norm) are the same solution.
  Scalar duplicate_tol = Scalar(1e-6);
  std::uint32_t seed = 12345;
};

template <typename Scalar>
struct MonodromyResult
{
  // Base parameters p0; every solution satisfies F(x, p0) = 0.
  std::vector<std::complex<Scalar>> parameters;
  std::vector<std::vector<std::complex<Scalar>>> solutions;
  int loops;
  int paths_tracked;
};

// Solves a parametric family F(x, p) = 0 at a generic parameter point by
// monodromy: starting from one known pair (x0, p0) with F(x0, p0) = 0, every
// known solution is tracked around triangular loops p0 -> p1 -> p2 -> p0 with
// random complex p1, p2. Since the solutions over p0 are permuted by the
// monodromy group, the loop endpoints are again solutions over p0, and (for
// an irreducible family) eventually every one of them is reached. Loops are
// kept and every new solution is tracked around all of them, so each round
// computes the orbit of x0 under the loops generated so far; a new random
// loop is added per round.
//
// p0 should be generic (e.g. random complex): for a start pair, it is
// usually easiest to pick x0 at random and solve F(x0, p) = 0 for p, which
// is linear if F is linear in its parameters.
template <typename Scalar>
class MonodromySolver
{
public:
  using Complex = std::complex<Scalar>;
  using Point = std::vector<Complex>;

  explicit MonodromySolver(MonodromyOptions<Scalar> options = {}) : options_(options) {}

  MonodromyResult<Scalar> solve(const System<Scalar>& F,
                                const std::vector<Variable>& parameters,
                                const Point& x0, const Point& p0) const
  {
    std::mt19937 rng(options_.seed);
    ParameterHomotopy<Scalar> H(F, parameters, p0, p0);
    EulerPredictor<Scalar> predictor;
    Tracker<Scalar> tracker(H, predictor, options_.corrector_tol, options_.corrector_max_iters,
                            options_.dt_init, options_.dt_min, options_.dt_max,
                            options_.max_steps, options_.max_norm);

    MonodromyResult<Scalar> result{p0, {x0}, 0, 0};
    std::vector<std::vector<Point>> loops;
    // tracked[i] = how many of `loops` solutions[i] has been tracked around.
    std::vector<size_t> tracked = {0};
    int loops_no_progress = 0;

    while (!done(result) && result.loops < options_.max_loops &&
           loops_no_progress < options_.max_loops_no_progress)
    {
      loops.push_back({p0, random_parameters(p0.size(), rng),
                       random_parameters(p0.size(), rng), p0});
      ++result.loops;
      size_t known = result.solutions.size();

      // Track every solution, including ones found along the way, around
      // every loop so far. Afterwards the solution set is closed under all
      // loops, so a round without growth means the new loop added nothing.
      for (size_t i = 0; i < result.solutions.size() && !done(result); ++i)
      {
        while (tracked[i] < loops.size() && !done(result))
        {
          Point x = result.solutions[i];
          ++result.paths_tracked;
          if (track_loop(x, loops[tracked[i]++], H, tracker) &&
              !contains(result.solutions, x))
          {
            result.solutions.push_back(std::move(x));
            tracked.push_back(0);
          }
        }
      }

      loops_no_progress = result.solutions.size() > known ? 0 : loops_no_progress + 1;
    }

    return result;
  }

private:
  bool done(const MonodromyResult<Scalar>& result) const
  {
    return options_.target_solutions_count > 0 &&
           static_cast<int>(result.solutions.size()) >= options_.target_solutions_count;
  }

  static bool track_loop(Point& x, const std::vector<Point>& loop,
                         ParameterHomotopy<Scalar>& H, const Tracker<Scalar>& tracker)
  {
    for (size_t k = 0; k + 1 < loop.size(); ++k)
    {
      H.set_parameters(loop[k], loop[k + 1]);
      auto r = tracker.track(std::move(x), Scalar(0), Scalar(1));
      if (!r.success)
        return false;
      x = std::move(r.solution);
    }
    return true;
  }

  bool contains(const std::vector<Point>& solutions, const Point& x) const
  {
    return std::any_of(solutions.begin(), solutions.end(), [&](const Point& s) {
      Scalar dist(0);
      for (size_t j = 0; j < x.size(); ++j)
        dist = std::max(dist, std::abs(s[j] - x[j]));
      return dist < options_.duplicate_tol;
    });
  }

  static Point random_parameters(size_t n, std::mt19937& rng)
  {
    std::normal_distribution<Scalar> normal(Scalar(0), Scalar(1));
    Point p(n);
    for (auto& pk : p)
      pk = Complex(normal(rng), normal(rng));
    return p;
  }

  MonodromyOptions<Scalar> options_;
};

} // namespace hc
