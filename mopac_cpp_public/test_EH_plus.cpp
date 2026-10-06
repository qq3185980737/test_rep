// test_EH_plus.cpp — numeric check for the third-generation H-bond correction.
#include "H_bond_correction_EH_plus.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "funcon_C.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;

double reada(const std::string&, int) { return 0.0; }  // stub (reada.F90)

int main() {
    int ok = 1;
    // Geometry:  D(0,0,0) H(1,0,0) A(3,0,0) linear
    //   D neighbour  RA=4=(0,1,0), RB=RC=6=(-1,0,0)
    //   A neighbour  AA=5=(3,1,0), AB=AC=7=(3,0,-1)
    // hblist row: 1=D 2=RA 3=RB 4=RC 5=A 6=AA 7=AB 8=AC 9=H 10=energy
    numat = 7;
    nat.assign(8, 0); nat[1] = 8; nat[2] = 1; nat[3] = 8; nat[4] = 1; nat[5] = 1; nat[6] = 1; nat[7] = 1;
    coord.assign(4, std::vector<double>(8, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.0; coord[2][2] = 0.0; coord[3][2] = 0.0;
    coord[1][3] = 3.0; coord[2][3] = 0.0; coord[3][3] = 0.0;
    coord[1][4] = 0.0; coord[2][4] = 1.0; coord[3][4] = 0.0;
    coord[1][5] = 3.0; coord[2][5] = 1.0; coord[3][5] = 0.0;
    coord[1][6] = -1.0; coord[2][6] = 0.0; coord[3][6] = 0.0;
    coord[1][7] = 3.0; coord[2][7] = 0.0; coord[3][7] = -1.0;

    const int max_h_bonds = 4;
    int hblist[4 * 10] = {0};
    int nrbondsa[4] = {1, 1, 1, 1};
    int nrbondsb[4] = {1, 1, 1, 1};
    hblist[0*4+1] = 1; hblist[1*4+1] = 4; hblist[2*4+1] = 6; hblist[3*4+1] = 6;
    hblist[4*4+1] = 3; hblist[5*4+1] = 5; hblist[6*4+1] = 7; hblist[7*4+1] = 7;
    hblist[8*4+1] = 2; hblist[9*4+1] = 0;

    const double a02 = a0 * a0;
    const double h2k = ev * fpc_9;

    // ---- non-PM7 branch ----
    method_pm7 = false;
    double E = EH_plus(1, hblist, max_h_bonds, nrbondsa, nrbondsb);
    // scale_c = -0.12*a0^2 ; angle_cos=1, angle2_cos=cos(2pi/3-pi/2)=0.8660,
    // torsion_cos=1 (RB==RC), angle2_cos_new=0.8660, torsion_cos_new=1 (AB==AC)
    double xy = 3.0;              // D-A distance
    double xc = 1.0;              // min(HA, HD) = 1.0
    double dmp = 1.0 - 1.0 / (1.0 + std::exp(-60.0 * (xc / 1.2 - 1.0)));
    dmp = dmp / (1.0 + std::exp(-100.0 * (xy / 2.4 - 1.0)));
    dmp = dmp * (1.0 - 1.0 / (1.0 + std::exp(-10.0 * (xy / 7.0 - 1.0))));
    double a2c = std::cos(2.0 * pi / 3.0 - pi / 2.0);
    double want = (-0.12 * a02) / (xy * xy) * 1.0 * a2c * a2c * 1.0 * a2c * a2c * 1.0 * h2k * dmp;
    if (std::fabs(E - want) > 1e-6 * std::max(1.0, std::fabs(want))) {
        std::printf("FAIL non-PM7: got %g want %g\n", E, want); ok = 0;
    }
    if (E > 0) { std::printf("FAIL sign non-PM7\n"); ok = 0; }

    // ---- PM7 branch (includes O-O short-range term) ----
    method_pm7 = true;
    E = EH_plus(1, hblist, max_h_bonds, nrbondsa, nrbondsb);
    double sc = -0.098822 * a02;
    double xy2 = 3.0;
    double xc2 = 1.0;
    double dm2 = 1.0;
    // XY_dist = max(ha,hb) - xc = 2-1 = 1 > 0.5 -> damped form
    dm2 = 1.0 - 1.0 / (1.0 + std::exp(-60.0 * (xc2 / 1.2 - 1.0)));
    dm2 = dm2 / (1.0 + std::exp(-100.0 * (xy2 / 2.4 - 1.0)));
    dm2 = dm2 * (1.0 - 1.0 / (1.0 + std::exp(-10.0 * (xy2 / 7.0 - 1.0))));
    double w2 = sc / (xy2 * xy2) * 1.0 *
                (1.0 - std::pow(1.0 - a2c * 1.0 * a2c * 1.0, 2)) * h2k * dm2;
    double short_o = -2.5 * std::exp(-80.0 * std::pow(std::max(xy2 - 2.67, 0.0), 2)) * 1.0;
    w2 += short_o;
    if (std::fabs(E - w2) > 1e-6 * std::max(1.0, std::fabs(w2))) {
        std::printf("FAIL PM7: got %g want %g\n", E, w2); ok = 0;
    }

    // ---- empty slot marker -> 0 ----
    hblist[9*4+1] = -666;
    E = EH_plus(1, hblist, max_h_bonds, nrbondsa, nrbondsb);
    if (E != 0.0) { std::printf("FAIL marker: %g\n", E); ok = 0; }
    hblist[9*4+1] = 0;

    // ---- angle_cos <= 0 -> 0: bend H well off line (D-H-A angle < 90 deg) ----
    coord[1][2] = 2.9; coord[2][2] = 1.0;
    E = EH_plus(1, hblist, max_h_bonds, nrbondsa, nrbondsb);
    if (E != 0.0) { std::printf("FAIL angle_cos<=0: %g\n", E); ok = 0; }

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
