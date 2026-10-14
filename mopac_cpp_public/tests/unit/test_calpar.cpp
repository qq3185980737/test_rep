// test_calpar.cpp
#include <cmath>
#include <cstdio>

#include "calpar.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

int main() {
    using namespace parameters_C;
    funcon_C::ev = 27.211386;
    molkst_C::keywrd = " PM3 ";

    // Hydrogen (i=1): only am[1]=gss[1]/ev.
    gss[1] = 5.0;
    // Carbon (i=6): s2 p2.
    ios[6] = 2; iop[6] = 2; iod[6] = 0;
    zs[6] = 0.6; zp[6] = 0.5;
    gss[6] = 10.0; gsp[6] = 11.0; gpp[6] = 12.0; gp2[6] = 4.0; hsp[6] = 1.0;
    uss[6] = -20.0; upp[6] = -10.0; udd[6] = 0.0;

    calpar();

    bool r1 = std::fabs(am[1] - 5.0 / funcon_C::ev) < 1e-9;
    std::printf("H am=%.4f (expect %.4f) %s\n", am[1], 5.0 / funcon_C::ev,
                r1 ? "PASS" : "FAIL");

    bool finite = std::isfinite(am[6]) && std::isfinite(ad[6]) && std::isfinite(aq[6]);
    bool r2 = finite && am[6] == 10.0 / funcon_C::ev && ad[6] > 0 && aq[6] > 0;
    std::printf("C am=%.4f ad=%.4f aq=%.4f eisol=%.4f %s\n", am[6], ad[6], aq[6],
                eisol[6], r2 ? "PASS" : "FAIL");

    return (r1 && r2) ? 0 : 1;
}
