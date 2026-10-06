// dot.cpp — C++ translation of MOPAC 2016 "dot.F90".
// Fortran x(n), y(n) are 1-based; convention keeps padding at [0].
#include "dot.h"

double dot(const std::vector<double>& x, const std::vector<double>& y, int n) {
    double s = 0.0;
    for (int i = 1; i <= n; ++i) s += x[i] * y[i];
    return s;
}
