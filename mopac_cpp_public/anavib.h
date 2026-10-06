// anavib.h — C++ translation of MOPAC 2016 "anavib.F90".
#pragma once
#include <vector>

// anavib: analyse vibrations and print description of normal modes.
// Arrays keep Fortran 1-based indexing (index 0 is padding).
//   eigs(n3), dipt(n3), hess(packed), rij(packed), f(packed): double, 1-based.
//   vibs(n3, n3): normal-mode eigenvectors, column-major (vibs(row, mode)).
void anavib(std::vector<double>& eigs, std::vector<double>& dipt, int n3,
            std::vector<std::vector<double>>& vibs, std::vector<double>& rij,
            int nv, std::vector<double>& hess, std::vector<double>& f);
