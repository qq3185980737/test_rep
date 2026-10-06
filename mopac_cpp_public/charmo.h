// charmo.h — C++ translation of MOPAC 2016 "charmo.F90".
#pragma once
#include <vector>

// charmo: symmetry transform of a MO under operation ioper.
double charmo(const std::vector<std::vector<double>>& vects,
              const std::vector<int>& ntype, int jorb, int ioper,
              const std::vector<std::vector<double>>& r, int nvecs, bool& first);
