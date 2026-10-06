// dist2.cpp — C++ translation of MOPAC 2016 "dist2.F90".

#include "dist2.h"

double dist2(const std::array<double,3>& a, const std::array<double,3>& b) {
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return dx*dx + dy*dy + dz*dz;
}
double dot1(const std::array<double,3>& a, const std::array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
