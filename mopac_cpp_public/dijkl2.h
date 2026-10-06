// dijkl2.h — C++ translation of MOPAC 2016 "dijkl2.F90".
#pragma once
#include <vector>

// Relax 2e integrals in MO basis into xy(4D).
void dijkl2(const std::vector<std::vector<double>>& dc);
