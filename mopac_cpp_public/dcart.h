// dcart.h — C++ translation of MOPAC 2016 "dcart.F90".
#pragma once
#include <vector>

// Smoothly-truncated -1/r^2 derivative used by point-charge term.
double derp(double r);

// Cartesian derivatives by finite differences (skeleton; dhc/delsta external).
void dcart(std::vector<std::vector<double>>& coord,
           std::vector<std::vector<double>>& dxyz);
