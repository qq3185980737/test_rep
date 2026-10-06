// test_bangle.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "bangle.h"
#include "molkst_C.h"

int main() {
    molkst_C::id = 0;
    // 3 atoms: i=1 at (1,0,0), j=2 at (0,0,0), k=3 at (0,1,0) -> 90 deg
    std::vector<std::vector<double>> xyz(4, std::vector<double>(4, 0.0));
    xyz[1][1] = 1.0; xyz[1][2] = 0.0; xyz[1][3] = 0.0;
    xyz[2][1] = 0.0; xyz[2][2] = 0.0; xyz[2][3] = 1.0;

    double ang;
    bangle(xyz, 1, 2, 3, ang);
    bool r1 = std::fabs(ang - 1.57079632679) < 1e-9;
    std::printf("angle i-j-k = %.6f rad (expect pi/2) %s\n", ang, r1 ? "PASS" : "FAIL");

    // straight line: k at (-1,0,0) -> 180 deg
    xyz[1][3] = -1.0; xyz[2][3] = 0.0;
    bangle(xyz, 1, 2, 3, ang);
    bool r2 = std::fabs(ang - 3.14159265359) < 1e-9;
    std::printf("angle straight = %.6f rad (expect pi) %s\n", ang, r2 ? "PASS" : "FAIL");

    return (r1 && r2) ? 0 : 1;
}
