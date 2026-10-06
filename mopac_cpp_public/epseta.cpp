// epseta.cpp — C++ translation of MOPAC 2016 "epseta.F90".

#include "epseta.h"

#include <algorithm>
#include <cfloat>
#include <cmath>

void epseta(double& eps, double& eta) {
    eps = std::max(DBL_MIN, 1e-39);
    eta = std::max(DBL_EPSILON, 1e-17);
}
