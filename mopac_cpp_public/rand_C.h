// rand_C.h — C++ translation of MOPAC 2016 "rand_C.F90": PCG64 RNG +
// Box-Muller normal transform (as added to this source tree).
#pragma once

#include <cstdint>
#include <vector>

namespace rand_C {
void init_random(int seed = 0);
std::uint64_t pcg64_next();
double pcg64_next_double();
void coupled_lorenz_generate(double x0, double y0, double z0, int T,
                             std::vector<double>& y);
void chaos_normal(double& z, const std::vector<double>& urand, int& ptr);
}  // namespace rand_C
