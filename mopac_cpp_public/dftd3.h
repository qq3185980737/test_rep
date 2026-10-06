// dftd3.h — C++ translation of MOPAC 2016 "dftd3.F90".
#pragma once
#include <vector>

// Grimme D3 dispersion correction. Returns dispersion energy (kcal/mol).
double dftd3(bool l_grad, std::vector<std::vector<double>>& dxyz);
