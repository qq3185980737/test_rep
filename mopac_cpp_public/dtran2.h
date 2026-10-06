// dtran2.h — C++ translation of MOPAC 2016 "dtran2.F90".
#pragma once
#include <vector>

// Transform d-orbital rotation: r(3,3) rotation -> t(5,5,ioper).
void dtran2(std::vector<std::vector<double>>& r,
            std::vector<std::vector<std::vector<double>>>& t, int ioper);
