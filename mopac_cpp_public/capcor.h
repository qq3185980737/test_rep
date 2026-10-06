// capcor.h — C++ translation of MOPAC 2016 "capcor.F90".
#pragma once
#include <vector>

// capcor: correction to electronic energy due to capped bonds (atom type 102).
double capcor(const std::vector<int>& nat, const std::vector<int>& nfirst,
              const std::vector<int>& nlast, const std::vector<double>& p,
              const std::vector<double>& h);
