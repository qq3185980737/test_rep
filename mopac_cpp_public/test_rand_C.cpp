// test_rand_C.cpp — rand_C: deterministic PCG64 sequence + Box-Muller stats.
#include "rand_C.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace rand_C;

int main() {
    bool ok = true;
    // Deterministic seed: sequence must be reproducible.
    init_random(12345);
    double a1 = pcg64_next_double(), b1 = pcg64_next_double();
    init_random(12345);
    double a2 = pcg64_next_double(), b2 = pcg64_next_double();
    if (a1 != a2 || b1 != b2) { std::printf("FAIL reproducibility\n"); ok = false; }
    if (!(a1 >= 0.0 && a1 < 1.0)) { std::printf("FAIL range a1=%g\n", a1); ok = false; }
    // Different seeds give different streams
    init_random(54321);
    double c1 = pcg64_next_double();
    if (c1 == a1) { std::printf("FAIL seed diversity\n"); ok = false; }
    // coupled_lorenz_generate fills T entries 1..T
    std::vector<double> y;
    coupled_lorenz_generate(0, 0, 0, 5, y);
    if (y.size() != 6) { std::printf("FAIL y size\n"); ok = false; }
    for (int i = 1; i <= 5; ++i) if (!(y[i] >= 0.0 && y[i] < 1.0)) ok = false;
    // chaos_normal Box-Muller: sample mean ~0, std ~1
    std::vector<double> u;
    coupled_lorenz_generate(0, 0, 0, 10002, u);
    int ptr = 1;
    double sum = 0, sumsq = 0;
    const int N = 5000;
    for (int i = 0; i < N; ++i) {
        double z;
        chaos_normal(z, u, ptr);
        sum += z;
        sumsq += z * z;
    }
    double mean = sum / N, var = sumsq / N - mean * mean;
    if (std::fabs(mean) > 0.1 || std::fabs(var - 1.0) > 0.15) {
        std::printf("FAIL normal stats mean=%g var=%g\n", mean, var);
        ok = false;
    }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
