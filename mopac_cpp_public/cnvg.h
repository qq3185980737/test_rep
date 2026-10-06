// cnvg.h — C++ translation of MOPAC 2016 "cnvg.F90".
#pragma once
#include <vector>

// Two-point interpolation for density-matrix convergence.
void cnvg(std::vector<double>& pnew, std::vector<double>& p,
          std::vector<double>& p1, int niter, double& pl);
