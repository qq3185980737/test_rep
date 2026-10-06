// deriv.h — C++ translation of MOPAC 2016 "deriv.F90".
#pragma once
#include <vector>

// Compute gradient wrt internal coordinates (skeleton; external routines stubbed).
void deriv(const std::vector<std::vector<double>>& geo,
           std::vector<double>& gradnt);
