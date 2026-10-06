// diag.h — C++ translation of MOPAC 2016 "diag.F90".
#pragma once
#include <vector>

// Fast pseudo-diagonalization: rotate old eigenvectors to block-diagonalize.
void diag(const std::vector<double>& fao,
          std::vector<std::vector<double>>& vector,
          int nocc, const std::vector<double>& eig, int mdim, int n);
