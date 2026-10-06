// dipole.h — C++ translation of MOPAC 2016 "dipole.F90".
#pragma once
#include <vector>

// Dipole moment (Debye). p = packed density, coord(3,n).
double dipole(const std::vector<double>& p,
              std::vector<std::vector<double>>& coord,
              std::vector<double>& dipvec, int mode);
