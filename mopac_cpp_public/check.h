// check.h — C++ translation of MOPAC 2016 "check.F90".
#pragma once
#include <vector>

// check: renormalize localized molecular orbitals.
void check(int nvec, const std::vector<int>& nnc, const std::vector<int>& nc,
           const std::vector<int>& icvec, const std::vector<int>& iorbs,
           const std::vector<int>& ncvec, std::vector<double>& cvec);
