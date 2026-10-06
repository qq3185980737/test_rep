// exchng.cpp — C++ translation of MOPAC 2016 "exchng.F90".

#include "exchng.h"

void exchng(double a, double& b, double c, double& d,
            double t, double& q, const std::vector<double>& x,
            std::vector<double>& y, int n) {
    b = a;
    d = c;
    q = t;
    for (int i = 1; i <= n; ++i) y[i] = x[i];
}
