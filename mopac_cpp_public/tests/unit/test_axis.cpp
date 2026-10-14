// test_axis.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "axis.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "to_screen_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace funcon_C;

int main() {
    // H2 along x: atoms at +/-0.37 A (r=0.74 A), mass 1.008 amu.
    numat = 2;
    mol_weight = 2.016;
    numcal = 1;
    atmass.assign(3, 0.0); atmass[1] = 1.008; atmass[2] = 1.008;
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[1][1] = -0.37; coord[1][2] = 0.37;

    double a, b, c;
    double evec[4][4] = {};
    axis(a, b, c, evec);

    // Expected: Ixx=0, Iyy=Izz = 2*1.008*0.37^2 = 0.2760 amu-A^2.
    double Izz = 2.0 * 1.008 * 0.37 * 0.37;
    double const2 = funcon_C::fpc_6 * funcon_C::fpc_10 * 1.0e16 /
                    (8.0 * funcon_C::pi * funcon_C::pi * funcon_C::fpc_8);
    double rot_expected = const2 / Izz;
    std::printf("a=%.4f b=%.4f c=%.4f  (expect ~%.3f %.3f 0)\n",
                a, b, c, rot_expected, rot_expected);
    bool r1 = std::fabs(a - rot_expected) < 0.05 && std::fabs(b - rot_expected) < 0.05 &&
              std::fabs(c) < 1.0;
    // evec orthonormal
    bool r2 = true;
    for (int i = 1; i <= 3 && r2; ++i) {
        double nrm = evec[1][i]*evec[1][i] + evec[2][i]*evec[2][i] + evec[3][i]*evec[3][i];
        if (std::fabs(nrm - 1.0) > 1e-9) r2 = false;
    }
    std::printf("evec orthonormal %s\n", r2 ? "PASS" : "FAIL");
    std::printf("rot constants %s\n", r1 ? "PASS" : "FAIL");
    return (r1 && r2) ? 0 : 1;
}
