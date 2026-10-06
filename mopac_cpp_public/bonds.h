// bonds.h — C++ translation of MOPAC 2016 "bonds.F90".
#pragma once
#include <vector>

// bonds: compute and print bond orders and valencies.
void bonds();

// dopen: ROHF open-shell density matrix from eigenvector matrix c(mdim,mdim).
// sdm(lower triangle, packed) = frac * sum_{nl1..nu1} c(i,n)*c(j,n).
void dopen(const std::vector<std::vector<double>>& eigvec, int mdim, int nor,
           int ndubl, int nsingl, double frac_in, std::vector<double>& sdm);
