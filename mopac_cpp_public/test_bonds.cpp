// test_bonds.cpp
#include <cstdio>
#include <vector>

#include "bonds.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

int main() {
    // c(i,n) = i*n for a 2x2 matrix. Active MO n=2 (ndubl=1, nsingl=2).
    int norbs = 2;
    std::vector<std::vector<double>> c(norbs + 1, std::vector<double>(norbs + 1, 0.0));
    for (int i = 1; i <= norbs; ++i)
        for (int n = 1; n <= norbs; ++n) c[i][n] = i * n;

    int ndubl = 1, nsingl = 2;  // nl1=2, nu1=2
    double fract = 1.0;
    std::vector<double> sdm((norbs * (norbs + 1)) / 2 + 1, 0.0);
    dopen(c, norbs, norbs, ndubl, nsingl, fract, sdm);

    // Expected: sdm(1,1)=c(1,2)^2=4; (2,1)=c(2,2)*c(1,2)=8; (2,2)=c(2,2)^2=16
    bool r1 = sdm[1] == 4.0 && sdm[2] == 8.0 && sdm[3] == 16.0;
    std::printf("sdm=(%.1f,%.1f,%.1f) %s\n", sdm[1], sdm[2], sdm[3], r1 ? "PASS" : "FAIL");

    // With fract=0.5, all values halve.
    dopen(c, norbs, norbs, ndubl, nsingl, 0.5, sdm);
    bool r2 = sdm[1] == 2.0 && sdm[3] == 8.0;
    std::printf("fract=0.5 sdm=(%.1f,%.1f) %s\n", sdm[1], sdm[3], r2 ? "PASS" : "FAIL");

    return (r1 && r2) ? 0 : 1;
}
