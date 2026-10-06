// diat.h — C++ translation of MOPAC 2016 "diat.F90".
#pragma once
#include <vector>

// Generalized Slater-type orbital overlap (pure); bfn external.
double ss_overlap(int na, int nb, int la, int lb, int m,
                  double ua, double ub, double r1, double a0);

// Diatomic overlap driver (skeleton).
void diat(int ni, int nj, const std::vector<double>& xj,
          std::vector<std::vector<double>>& di);
