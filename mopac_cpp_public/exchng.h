// exchng.h — C++ translation of MOPAC 2016 "exchng.F90".
#pragma once
#include <vector>
void exchng(double a, double& b, double c, double& d,
            double t, double& q, const std::vector<double>& x,
            std::vector<double>& y, int n);
