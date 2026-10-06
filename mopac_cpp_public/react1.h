// react1.h — C++ translation of MOPAC 2016 "react1.F90".
#pragma once
#include <vector>
void react1();
void dock(std::vector<std::vector<double>>& geoa,
          std::vector<std::vector<double>>& geo, double& dist);
