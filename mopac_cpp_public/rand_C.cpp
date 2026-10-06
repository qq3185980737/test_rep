// rand_C.cpp — C++ translation of MOPAC 2016 "rand_C.F90".
// PCG64-style generator (32-bit output) + Box-Muller normal transform.
#include "rand_C.h"

#include <cmath>
#include <cstdio>

namespace rand_C {

namespace {
std::uint64_t state = 0;
std::uint64_t inc = 0;
constexpr std::uint64_t MULTIPLIER = 6364136223846793005ULL;
bool initialized = false;
constexpr double TWO53 = 9007199254740992.0;  // 2^53
}  // namespace

void init_random(int seed) {
    std::uint64_t seed_val = static_cast<std::uint64_t>(seed);
    inc = 2ULL + 1ULL;  // must be odd (Fortran: rng%inc = (1_8 * 2_8) + 1_8)
    state = 0;
    (void)pcg64_next();              // advance once
    state += seed_val;
    (void)pcg64_next();              // advance once more
    initialized = true;
}

std::uint64_t pcg64_next() {
    const std::uint64_t oldstate = state;
    state = oldstate * MULTIPLIER + inc;
    // PCG-XSH-RR with 32-bit output
    const std::uint32_t xorshifted =
        static_cast<std::uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
    const std::uint32_t rot = static_cast<std::uint32_t>(oldstate >> 59u);
    std::uint32_t value;
    if (rot == 0) {
        value = xorshifted;
    } else {
        value = (xorshifted >> rot) | (xorshifted << (32u - rot));
    }
    return static_cast<std::uint64_t>(value);
}

double pcg64_next_double() {
    // pcg64_next() yields 32 bits; draw twice and splice to a 53-bit mantissa
    // (the Fortran comment intends "take the high 53 bits" of the stream).
    const std::uint64_t hi = pcg64_next();
    const std::uint64_t lo = pcg64_next();
    const std::uint64_t bits = (hi << 21u) | (lo >> 11u);
    double value = static_cast<double>(bits) / TWO53;
    if (value >= 1.0) value = 1.0 - 1.0e-15;
    if (value < 0.0) value = 0.0;
    return value;
}

void coupled_lorenz_generate(double /*x0*/, double /*y0*/, double /*z0*/,
                             int T, std::vector<double>& y) {
    if (!initialized) init_random();
    y.assign(static_cast<size_t>(T) + 1, 0.0);  // 1-based
    for (int i = 1; i <= T; ++i) y[i] = pcg64_next_double();
}

void chaos_normal(double& z, const std::vector<double>& urand, int& ptr) {
    const double pi = 4.0 * std::atan(1.0);
    if (ptr + 1 > static_cast<int>(urand.size()) - 1) {
        std::printf("错误: chaos_normal 指针越界, ptr=%d size=%d\n", ptr,
                    static_cast<int>(urand.size()));
        z = 0.0;
        return;
    }
    double u1 = urand[ptr];
    double u2 = urand[ptr + 1];
    ptr += 2;
    u1 = std::max(u1, 1.0e-12);
    u1 = std::min(u1, 1.0 - 1.0e-12);
    u2 = std::max(u2, 1.0e-12);
    u2 = std::min(u2, 1.0 - 1.0e-12);
    z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * pi * u2);
}

}  // namespace rand_C
