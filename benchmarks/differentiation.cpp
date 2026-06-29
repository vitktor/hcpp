#include <benchmark/benchmark.h>
#include <hc/core/differentiation.hpp>
#include <hc/core/system.hpp>

using namespace hc;

static void DifferentiateUnivariate(benchmark::State& state) {
    Variable x("x");
    Polynomial<double> p({1.0, 2.0, 1.0, 1.0}, {{3}, {2}, {1}, {0}}, {x});
    for (auto _ : state) {
        benchmark::DoNotOptimize(differentiate(p, x));
    }
}
BENCHMARK(DifferentiateUnivariate);

static void DifferentiateMultivariate(benchmark::State& state) {
    Variable x("x"), y("y");
    Polynomial<double> p({1.0, 1.0, 1.0}, {{2, 0}, {1, 1}, {0, 2}}, {x, y});
    for (auto _ : state) {
        benchmark::DoNotOptimize(differentiate(p, x));
    }
}
BENCHMARK(DifferentiateMultivariate);

static void Jacobian2x2(benchmark::State& state) {
    Variable x("x"), y("y");
    Polynomial<double> f1({1.0, 1.0, -1.0}, {{2, 0}, {0, 2}, {0, 0}}, {x, y});
    Polynomial<double> f2({1.0, -1.0}, {{1, 0}, {0, 1}}, {x, y});
    for (auto _ : state) {
        benchmark::DoNotOptimize(System<double>({f1, f2}, {x, y}));
    }
}
BENCHMARK(Jacobian2x2);

BENCHMARK_MAIN();
