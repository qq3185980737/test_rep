// chrge.h — C++ translation of MOPAC 2016 "chrge.F90".
#pragma once
#include <vector>

// Compute atom electron densities from packed density matrix p.
void chrge(const std::vector<double>& pin, std::vector<double>& q);
