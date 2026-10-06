// dihed.h — C++ translation of MOPAC 2016 "dihed.F90".
#pragma once
#include <vector>

// Dihedral angle (radians) between atoms i,j,k,l in xyz(3,n).
void dihed(const std::vector<std::vector<double>>& xyz,
           int i, int j, int k, int l, double& angle);
