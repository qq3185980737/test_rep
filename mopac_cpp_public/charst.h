// charst.h — C++ translation of MOPAC 2016 "charst.F90".
#pragma once
#include <vector>

// charst: character of CI state istate under operation ioper.
// istate < 0 resets cached microstate state and returns 0.
double charst(const std::vector<std::vector<double>>& vects,
              const std::vector<int>& ntype, int istate, int ioper,
              const std::vector<std::vector<double>>& r, int nvecs, bool& first);
