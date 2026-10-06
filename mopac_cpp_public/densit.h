// densit.h — C++ translation of MOPAC 2016 "densit.F90".
#pragma once
#include <vector>

// Packed lower-triangular density matrix p from square eigenvector matrix c.
void densit(const std::vector<std::vector<double>>& c, int mdim, int norbs,
            int nocc, double occ, int nfract, double fract,
            std::vector<double>& p, int mode);
