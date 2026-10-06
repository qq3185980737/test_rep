// deri1.h — C++ translation of MOPAC 2016 "deri1.F90".
#pragma once
#include <vector>

// Non-relaxed derivative wrt one Cartesian coordinate (skeleton).
void deri1(int number, double& grad,
           std::vector<double>& fv, int minear,
           std::vector<double>& fd,
           const std::vector<double>& scalar,
           std::vector<std::vector<double>>& work);
