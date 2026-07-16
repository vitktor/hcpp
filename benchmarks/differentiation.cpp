#include <benchmark/benchmark.h>
#include <hc/core/calculus.hpp>
#include <hc/core/system.hpp>

using namespace hc;

static void DifferentiateUnivariate(benchmark::State& state) {
    Variable x("x");
    auto p = pow(x, 3) + 2.0*pow(x, 2) + x + 1.0;
    for (auto _ : state) {
        benchmark::DoNotOptimize(differentiate(p, x));
    }
}
BENCHMARK(DifferentiateUnivariate);

static void DifferentiateMultivariate(benchmark::State& state) {
    Variable x("x"), y("y");
    auto p = pow(x, 2) + x*y + pow(y, 2);
    for (auto _ : state) {
        benchmark::DoNotOptimize(differentiate(p, x));
    }
}
BENCHMARK(DifferentiateMultivariate);

static void Jacobian2x2(benchmark::State& state) {
    Variable x("x"), y("y");
    auto f1 = pow(x, 2) + pow(y, 2) - 1.0;
    auto f2 = x - y;
    for (auto _ : state) {
        benchmark::DoNotOptimize(System<double>({f1, f2}, {x, y}));
    }
}
BENCHMARK(Jacobian2x2);

BENCHMARK_MAIN();
