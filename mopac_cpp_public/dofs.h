// dofs.h — C++ translation of MOPAC 2016 "dofs.F90".
#pragma once
#include <vector>

// Density-of-states histogram into dd(m). eref(mono3,n).
void dofs(std::vector<std::vector<double>>& eref, int mono3, int n,
          std::vector<double>& dd, int m, double bottom, double top);
