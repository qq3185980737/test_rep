// test_EC_plus_ER.cpp — numeric check for the PM6-DH2 H-bond correction.
#include "H_bond_correction_EC_plus_ER.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "funcon_C.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;

#include <string>
double reada(const std::string&, int) { return 0.0; }  // stub (reada.F90)

int main() {
    int ok = 1;
    // Water-acceptor geometry:
    //   D(0,0,0) - H(1,0,0) ... A(3,0,0)  linear -> angle_cos = 1
    //   A neighbours: 4=(3,1,0), 5=(3,cos1,sin1)  (both H -> water multiplier)
    numat = 5;
    nat.assign(6, 0); nat[1] = 8; nat[2] = 1; nat[3] = 8; nat[4] = 1; nat[5] = 1;
    nbonds.assign(6, 0); nbonds[1] = 1; nbonds[2] = 1; nbonds[3] = 2; nbonds[4] = 1; nbonds[5] = 1;
    ibonds.assign(16, std::vector<int>(6, 0));
    ibonds[1][1] = 2; ibonds[1][2] = 1;
    ibonds[1][3] = 4; ibonds[2][3] = 5;
    coord.assign(4, std::vector<double>(6, 0.0));
    coord[1][1] = 0.0; coord[2][1] = 0.0; coord[3][1] = 0.0;
    coord[1][2] = 1.0; coord[2][2] = 0.0; coord[3][2] = 0.0;
    coord[1][3] = 3.0; coord[2][3] = 0.0; coord[3][3] = 0.0;
    coord[1][4] = 3.0; coord[2][4] = 1.0; coord[3][4] = 0.0;
    coord[1][5] = 3.0; coord[2][5] = std::cos(1.0); coord[3][5] = std::sin(1.0);

    int d_list[8] = {0};
    int nd_list = 0;
    double EC = 0, ER = 0;

    // D=1 H=2 A=3, q1=q2=1
    double c = EC_plus_ER(1, 2, 3, 1.0, 1.0, EC, ER, d_list, nd_list);

    if (nd_list != 5) { std::printf("FAIL nd_list=%d\n", nd_list); ok = 0; }
    if (!(d_list[1] == 1 && d_list[2] == 2 && d_list[3] == 3 &&
          d_list[4] == 4 && d_list[5] == 5)) {
        std::printf("FAIL d_list: %d %d %d %d %d\n", d_list[1], d_list[2], d_list[3], d_list[4], d_list[5]);
        ok = 0;
    }

    // Hand-derived expectation for this geometry:
    //   angle_cos=1, angle2_cos=cos(pi*109.48/180 - pi/2)=0.942647,
    //   torsion_cos=cos(pi*54.74/180 - pi/2)=0.816612 (dihed=4.712389 -> pi/2),
    //   r = truncation(2,1.8,0.05)=2 -> r/a0=3.779448
    //   attraction = 0.76/ r^3 * fpc_9*ev, repulsion = 0.65*5^-r * fpc_9*ev
    double ra0 = 2.0 / a0;
    double unit = fpc_9 * ev;
    double attr = 0.76 / std::pow(ra0, 3) * unit;       // water multiplier
    double rep = 0.65 * std::pow(5.0, -ra0) * unit;
    double ac = 1.0, a2c = std::cos(pi * 109.48 / 180.0 - pi / 2.0);
    double tc = std::cos(pi * 54.74 / 180.0 - pi / 2.0);
    double e_EC = attr * ac * a2c * tc;
    double e_ER = rep * ac * a2c * tc;
    double e_c  = (attr + rep) * ac * a2c * tc;

    if (std::fabs(EC - e_EC) > 1e-6 * std::max(1.0, std::fabs(e_EC))) {
        std::printf("FAIL EC: got %g want %g\n", EC, e_EC); ok = 0;
    }
    if (std::fabs(ER - e_ER) > 1e-6 * std::max(1.0, std::fabs(e_ER))) {
        std::printf("FAIL ER: got %g want %g\n", ER, e_ER); ok = 0;
    }
    if (std::fabs(c - e_c) > 1e-6 * std::max(1.0, std::fabs(e_c))) {
        std::printf("FAIL c: got %g want %g\n", c, e_c); ok = 0;
    }
    if (EC < 0 || ER < 0) { std::printf("FAIL sign\n"); ok = 0; }

    // Sanity: q2<0 flips attraction only (repulsion is charge-independent)
    double EC2 = 0, ER2 = 0;
    int nd2 = 0;
    int dl2[8] = {0};
    double c2 = EC_plus_ER(1, 2, 3, 1.0, -1.0, EC2, ER2, dl2, nd2);
    double e_EC2 = -attr * ac * a2c * tc;
    double e_c2  = (-attr + rep) * ac * a2c * tc;
    if (std::fabs(EC2 - e_EC2) > 1e-6) { std::printf("FAIL EC sign flip: %g want %g\n", EC2, e_EC2); ok = 0; }
    if (std::fabs(c2 - e_c2) > 1e-6)   { std::printf("FAIL c sign flip: %g want %g\n", c2, e_c2); ok = 0; }

    // Degenerate: angle_cos < 0 (D-H-A angle < 90 deg) -> 0
    // H well off the D-A line and near A: D-H-A angle ~76 deg -> -cos<0
    coord[1][2] = 2.9; coord[2][2] = 1.0;
    double EC3 = 0, ER3 = 0;
    int nd3 = 0; int dl3[8] = {0};
    double c3 = EC_plus_ER(1, 2, 3, 1.0, 1.0, EC3, ER3, dl3, nd3);
    if (c3 != 0.0 || EC3 != 0.0 || ER3 != 0.0) { std::printf("FAIL angle_cos<0: %g %g %g\n", c3, EC3, ER3); ok = 0; }

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
