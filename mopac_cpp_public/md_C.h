// md_C.h — C++ translation of MOPAC 2016 "md_C.F90" (molecular dynamics globals).
#pragma once

#include <vector>

namespace md_C {
// Defaults (Fortran parameter constants)
constexpr int DEF_MDSTEP = 1000;
constexpr double DEF_DTMD = 0.2;
constexpr double DEF_TTARGET = 16000.0;
constexpr double KB = 3.166811563e-6;

// MD global arrays (3 x numat, 1-based)
extern std::vector<std::vector<double>> md_vel;
extern std::vector<std::vector<double>> md_force;
extern std::vector<std::vector<double>> md_vel_half;

// MD controls
extern int md_step_max;
extern double dt_md;
extern double T_target;
extern bool l_mdzero;
extern bool md_stop;
extern bool md_test;
}  // namespace md_C
