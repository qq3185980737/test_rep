// forsav.h — C++ translation of MOPAC 2016 "forsav.F90".
#pragma once
#include <vector>
void forsav(double& time, std::vector<std::vector<double>>& deldip,
            int& ipt, std::vector<double>& fmatrx, std::vector<double>& coord,
            int nvar, double& refh, std::vector<double>& evecs,
            int& jstart, std::vector<double>& fconst);
