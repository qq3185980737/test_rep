// fock2z.h — C++ translation of MOPAC 2016 "fock2z.F90".
#pragma once
#include <vector>

void focd2z(int iab, int jba,
            std::vector<double>& fii, std::vector<double>& fjj,
            std::vector<double>& fij,
            const std::vector<double>& pii, const std::vector<double>& pjj,
            const std::vector<double>& pij,
            const std::vector<double>& wj, const std::vector<double>& wk,
            bool diagonal, int& kr);
void fock2z(std::vector<double>& f, std::vector<double>& q,
            std::vector<double>& qe, const std::vector<double>& wj,
            const std::vector<double>& wk,
            std::vector<std::vector<double>>& ptot2, int mode, int ione);
