// md_C.cpp — C++ translation of MOPAC 2016 "md_C.F90" (molecular dynamics globals).
#include "md_C.h"

namespace md_C {
std::vector<std::vector<double>> md_vel;
std::vector<std::vector<double>> md_force;
std::vector<std::vector<double>> md_vel_half;

int md_step_max = DEF_MDSTEP;
double dt_md = DEF_DTMD;
double T_target = DEF_TTARGET;
bool l_mdzero = false;
bool md_stop = false;
bool md_test = true;
}  // namespace md_C
