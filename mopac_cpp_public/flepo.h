// flepo.h — C++ translation of MOPAC 2016 "flepo.F90".
#pragma once
#include <vector>

// Inverse-Hessian update (packed lower triangle, length nvar*(nvar+1)/2).
// dfp=true: DFP; false: BFGS. hesinv updated in place.
void flepo_update_H(int nvar, const std::vector<double>& xvar,
                    const std::vector<double>& gvar,
                    const std::vector<double>& gg,
                    double sy, double yhy, bool dfp,
                    std::vector<double>& hesinv);
void flepo(std::vector<double>& xparam, int nvar, double& funct1);
