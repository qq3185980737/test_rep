// dipole.cpp — C++ translation of MOPAC 2016 "dipole.F90".
// DIPOLE CALCULATES DIPOLE MOMENTS (ZDO: point-charge + one-center
// hybridization terms).
#include "dipole.h"

#include <cmath>
#include <string>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "to_screen_C.h"
#include "to_screen.h"

using namespace common_arrays_C;   // incl. nfirst, nlast, nat, atmass, q
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;
using namespace to_screen_C;
using chanel_C::iw;

double dipole(const std::vector<double>& pv,
              std::vector<std::vector<double>>& crd,
              std::vector<double>& dipvec, int mode) {
    static int icalcn = 0;
    static bool first = true, chargd = false, force = false;
    std::vector<double> center(4, 0.0);
    double sum, hyfsp, xt, hyfpd, dx, dy, dz, factdm;
    int i, j, ni, ia, ib, l, k, ll;
    first = (icalcn != numcal);
    icalcn = numcal;
    if (first) {
        sum = 0.0;
        for (i = 1; i <= numat; ++i) sum += q[i];
        chargd = std::fabs(sum) > 0.5;
        force = (keywrd.find("FORCE") != std::string::npos ||
                 keywrd.find(" THERMO") != std::string::npos ||
                 keywrd.find("IRC") != std::string::npos);
    }
    if (!force && chargd) {
        // RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE ORIGIN
        center.assign(4, 0.0);
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j)
                center[i] += atmass[j] * crd[i - 1][j];
        for (i = 1; i <= 3; ++i) center[i] /= mol_weight;
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j) crd[i - 1][j] -= center[i];
    }
    for (j = 0; j < 4; ++j)
        for (k = 0; k < 3; ++k) dip[j][k] = 0.0;
    for (i = 1; i <= numat; ++i) {
        ni = nat[i];
        ia = nfirst[i];
        ib = nlast[i];
        l = ib - ia;
        if (l > 0) {
            hyfsp = 2.0 * dd[ni] * a0 * fpc_8 * fpc_1 * 1.0e-10;
            for (j = 1; j <= 3; ++j) {
                k = ((ia + j) * (ia + j - 1)) / 2 + ia;
                dip[j - 1][1] -= hyfsp * pv[k];
            }
            // PD ONE-CENTER TERM
            if (ib - ia == 8) {
                xt = 1.0 / std::sqrt(3.0);
                hyfpd = 2.0 * ddp[5][ni] * a0 * fpc_8 * fpc_1 * 1.0e-10;
                ll = (ia + 5) * (ia + 4) / 2 + ia + 3;
                dx = pv[ll];
                ll = (ia + 4) * (ia + 3) / 2 + ia + 1;
                dx += pv[ll];
                ll = (ia + 8) * (ia + 7) / 2 + ia + 2;
                dx += pv[ll];
                ll = (ia + 6) * (ia + 5) / 2 + ia + 1;
                dx -= xt * pv[ll];
                ll = (ia + 7) * (ia + 6) / 2 + ia + 3;
                dy = pv[ll];
                ll = (ia + 4) * (ia + 3) / 2 + ia + 2;
                dy -= pv[ll];
                ll = (ia + 8) * (ia + 7) / 2 + ia + 1;
                dy += pv[ll];
                ll = (ia + 6) * (ia + 5) / 2 + ia + 2;
                dy -= xt * pv[ll];
                ll = (ia + 5) * (ia + 4) / 2 + ia + 1;
                dz = pv[ll];
                ll = (ia + 7) * (ia + 6) / 2 + ia + 2;
                dz += pv[ll];
                ll = (ia + 6) * (ia + 5) / 2 + ia + 3;
                dz += 2.0 * xt * pv[ll];
                dip[0][1] -= dx * hyfpd;
                dip[1][1] -= dy * hyfpd;
                dip[2][1] -= dz * hyfpd;
            }
        }
        // FPC(8)=SPEED OF LIGHT, FPC(1)=CHARGE ON ELECTRON
        factdm = fpc_8 * fpc_1 * 1.0e-10;
        for (j = 1; j <= 3; ++j)
            dip[j - 1][0] += q[i] * crd[j - 1][i] * factdm;
    }
    for (j = 1; j <= 3; ++j) dip[j - 1][2] = dip[j - 1][1] + dip[j - 1][0];
    for (j = 1; j <= 3; ++j)
        dip[3][j - 1] = std::sqrt(dip[0][j - 1] * dip[0][j - 1] +
                                  dip[1][j - 1] * dip[1][j - 1] +
                                  dip[2][j - 1] * dip[2][j - 1]);
    fprintf(stderr, "[DIPDBG] mode=%d q:", mode);
    for (i = 1; i <= numat; ++i) fprintf(stderr, " %+.6f", q[i]);
    fprintf(stderr, "\n[DIPDBG] crd:");
    for (i = 1; i <= numat; ++i)
        fprintf(stderr, " (%+.6f,%+.6f,%+.6f)", crd[0][i], crd[1][i], crd[2][i]);
    fprintf(stderr, "\n[DIPDBG] dip[0..3][0..2]: %.6f %.6f %.6f | %.6f %.6f %.6f | %.6f %.6f %.6f | %.6f %.6f %.6f\n",
        dip[0][0], dip[1][0], dip[2][0], dip[3][0],
        dip[0][1], dip[1][1], dip[2][1], dip[3][1],
        dip[0][2], dip[1][2], dip[2][2], dip[3][2]); fflush(stderr);
    if (force) {
        dipvec[1] = dip[0][2];
        dipvec[2] = dip[1][2];
        dipvec[3] = dip[2][2];
    }
    if (mode == 1) {
        if (iw > 0) {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                          " DIPOLE       X         Y         Z       TOTAL\n"
                          " POINT-CHG. %10.3f%10.3f%10.3f%10.3f\n"
                          " HYBRID   %10.3f%10.3f%10.3f%10.3f\n"
                          " SUM      %10.3f%10.3f%10.3f%10.3f",
                          dip[0][0], dip[1][0], dip[2][0], dip[3][0],
                          dip[0][1], dip[1][1], dip[2][1], dip[3][1],
                          dip[0][2], dip[1][2], dip[2][2], dip[3][2]);
            std::printf("%s\n", buf);   // DIPOLE table goes to the .out file
        }
    }
    // STORE DIPOLE MOMENT COMPONENTS IN UX,UY,UZ FOR USE IN ASSIGNING
    // CHARGES DETERMINED FROM THE ESP.
    ux = dip[0][2];
    uy = dip[1][2];
    uz = dip[2][2];
    if (!force && chargd) {
        // RESET ION'S POSITION AGAIN, SO THAT IT IS BACK WHERE IT STARTED
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j) crd[i - 1][j] += center[i];
    }
    return dip[3][2];
}
