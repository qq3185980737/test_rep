// bangle.h — C++ translation of MOPAC 2016 "bangle.F90".
#pragma once
#include <vector>

// bangle: angle (radians) at atom j between atoms i-j-k.
// xyz(3, atom) Cartesian coordinates (1-based).
void bangle(const std::vector<std::vector<double>>& xyz, int i, int j, int k,
            double& angle);
