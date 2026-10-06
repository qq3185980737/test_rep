// dijkl1.h — C++ translation of MOPAC 2016 "dijkl1.F90".
#pragma once
#include <vector>

// Two-electron integral derivatives over active MOs w.r.t. atom nati.
void dijkl1(const std::vector<std::vector<double>>& c, int n, int nati,
            const std::vector<double>& w,
            std::vector<double>& cij, std::vector<double>& wcij,
            std::vector<double>& ckl,
            std::vector<double>& xy);
