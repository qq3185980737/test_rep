// test_analyt.cpp
#pragma warning(disable: 4459)
#include <cmath>
#include <cstdio>
#include <vector>

#include "analyt.h"
#include "analyt_C.h"
#include "parameters_C.h"

using namespace analyt_C;
using namespace parameters_C;

int main() {
    // --- delri: H-H case (natorb <= 2 both) ---
    natorb[1] = 1; natorb[2] = 1;
    am[1] = 1.0; am[2] = 1.0;
    double dg[23] = {};
    double rr = 2.0, del1 = 2.0;  // along x
    delri(dg, 1, 2, rr, del1);
    double aee = 0.25 * (1.0 + 1.0) * (1.0 + 1.0);
    double ee = -rr / std::pow(std::sqrt(rr * rr + aee), 3);
    double a0 = 0.5291772083, ev = 27.2113834;
    double term = ev * del1 / (rr * a0 * a0);
    double expect = term * ee;
    bool ok = std::fabs(dg[1] - expect) < 1e-9;
    std::printf("delri HH dg[1]=%.6f expect=%.6f  %s\n", dg[1], expect, ok ? "PASS" : "FAIL");
    // dg(2..22) untouched in H-H early return
    ok &= (dg[2] == 0.0);

    // --- rotat: bond along +x, ix=1, idx=2 ---
    std::vector<std::vector<double>> coord(4, std::vector<double>(3, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;  // j=1
    coord[1][2] = 1.0; coord[2][2] = 0.0; coord[3][2] = 0.0;  // i=2
    rotat(coord, 2, 1, 1, 2, 1.0, 1.0);
    bool r1 = std::fabs(tx[1] - 1.0) < 1e-12 && std::fabs(tx[2]) < 1e-12 &&
              std::fabs(tx[3]) < 1e-12;
    std::printf("rotat tx=(%g,%g,%g) %s\n", tx[1], tx[2], tx[3], r1 ? "PASS" : "FAIL");
    // orthonormal: tx,ty,tz dot products ~0, norms ~1
    double dot_xy = tx[1]*ty[1]+tx[2]*ty[2]+tx[3]*ty[3];
    double nx = std::sqrt(tx[1]*tx[1]+tx[2]*tx[2]+tx[3]*tx[3]);
    bool r2 = std::fabs(dot_xy) < 1e-12 && std::fabs(nx - 1.0) < 1e-12;
    std::printf("rotat orthogonality %s\n", r2 ? "PASS" : "FAIL");

    return (ok && r1 && r2) ? 0 : 1;
}
