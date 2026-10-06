// compfg.h — C++ translation of MOPAC 2016 "compfg.F90".
// The SCF driver is a skeleton; xfac_value (core-core energy) is fully implemented.
#pragma once
#include <vector>

void compfg(const std::vector<double>& xparam, bool int_flag, double& escf,
            bool fulscf, std::vector<double>& grad, bool lgrad);

// Core-core repulsion used by xfac keyword.
double xfac_value();
