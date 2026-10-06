// deri0.h — C++ translation of MOPAC 2016 "deri0.F90".
#pragma once
#include <vector>

// Build diagonal dominant block and row scale factors for deri2 relaxation.
void deri0(const std::vector<double>& e, int n,
           std::vector<double>& scalar, std::vector<double>& diag,
           double fract, const std::vector<int>& nbo);
