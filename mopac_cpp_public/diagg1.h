// diagg1.h — C++ translation of MOPAC 2016 "diagg1.F90".
#pragma once
#include <vector>

// Build significant occupied-virtual LMO matrix elements fmo/ifmo (skeleton).
void diagg1(const std::vector<double>& fao, int nocc, int nvir,
            std::vector<double>& eigv, std::vector<double>& ws,
            std::vector<char>& latoms,
            std::vector<std::vector<int>>& ifmo, std::vector<double>& fmo,
            int fmo_dim, int& nij, int idiagg,
            std::vector<double>& avir, std::vector<double>& aocc,
            std::vector<double>& aov);
