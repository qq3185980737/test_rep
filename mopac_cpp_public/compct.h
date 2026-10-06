// compct.h — C++ translation of MOPAC 2016 "compct.F90".
#pragma once
#include <vector>

// Compress LMO data (ic / c / nc / iws) upward, dropping empty LMOs.
void compct(std::vector<int>& nncnew, std::vector<int>& ncnew,
            std::vector<int>& ncmnew, int itop, std::vector<int>& nc,
            std::vector<int>& ic, std::vector<int>& iws, int n01,
            std::vector<double>& c, int n02, int nmos, int idone,
            int& lb, int& mb, int icref, int inref);
