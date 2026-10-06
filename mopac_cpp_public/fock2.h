// fock2.h — C++ translation of MOPAC 2016 "fock2.F90".
#pragma once
#include <vector>

void fock1dorbs(std::vector<double>& f, const std::vector<double>& ptot,
                const std::vector<double>& pa,
                const std::vector<std::vector<double>>& w,
                int& kr, int ia, int ib, int ilim);
void fock2(std::vector<double>& f, const std::vector<double>& ptot,
           std::vector<double>& p, const std::vector<double>& w,
           const std::vector<double>& wj, const std::vector<double>& wk,
           int numat, const std::vector<int>& nfirst,
           const std::vector<int>& nlast, int mode);
