// dtrans.h — C++ translation of MOPAC 2016 "dtrans.F90".
#pragma once
#include <vector>

// Apply symmetry operation ioper to the 5 d-orbitals d.
void dtrans(std::vector<double>& d, int ioper, bool& first,
            const std::vector<std::vector<double>>& r);
