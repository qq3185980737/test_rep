// dipole_for_MOZYME.cpp — C++ translation of MOPAC 2016
// "dipole_for_MOZYME.F90".
#include "dipole_for_MOZYME.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "ijbo.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "to_screen.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;
using chanel_C::iw;

double dipole_for_MOZYME(std::vector<double>& dipvec, int mode) {
    static bool chargd = false, first = true, force = false;
    static int icalcn = 0;
    static double wtmol_loc = 0.0;
    static std::vector<double> hyf(108, 0.0);
    static std::vector<std::vector<double>> dip(4, std::vector<double>(3, 0.0));
    std::vector<double> center(4, 0.0);
    int i, j, k, l, ni;
    double sum, dx, dy, dz, hyfpd, xt;
    // CONST = c.C.a0.2
    const double fconst = fpc_8 * fpc_1 * a0 * 2.0 * 1e-10;
    first = (icalcn != numcal);
    icalcn = numcal;
    if (first) {
        hyf[1] = 0.0;
        for (i = 2; i <= 107; ++i) hyf[i] = fconst * dd[i];
        wtmol_loc = 0.0;
        sum = 0.0;
        for (i = 1; i <= numat; ++i) {
            wtmol_loc += ams[nat[i]];
            sum += q[i];
        }
        chargd = std::fabs(sum) > 0.5;
        force = (keywrd.find("FORCE") != std::string::npos ||
                 keywrd.find(" THERMO") != std::string::npos ||
                 keywrd.find("IRC") != std::string::npos);
    }
    if (chargd) {
        // RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE ORIGIN
        center.assign(4, 0.0);
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j)
                center[i] += ams[nat[j]] * coord[i-1][j];
        for (i = 1; i <= 3; ++i) center[i] /= wtmol_loc;
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j) coord[i-1][j] -= center[i];
    }
    for (j = 0; j < 4; ++j)
        for (k = 0; k < 3; ++k) dip[j][k] = 0.0;
    for (i = 1; i <= numat; ++i) {
        ni = nat[i];
        l = MOZYME_C::iorbs[i] - 1;
        if (l >= 3) {
            for (j = 1; j <= 3; ++j) {
                k = ijbo(i, i) + 1 + (j * (j + 1)) / 2;
                dip[j - 1][1] -= hyf[ni] * p[k];
            }
            if (l == 8) {
                xt = 1.0 / std::sqrt(3.0);
                hyfpd = 2.0 * ddp[5][ni] * a0 * fpc_8 * fpc_1 * 1e-10;
                k = ijbo(i, i);
                // x: <x|x|x2-y2> + <z|x|xz> - 1/root3<x|x|z2> + <y|x|xy>
                // y: -<y|y|x2-y2> - 1/root3<y|y|z2> + <z|y|yz> + <x|y|xy>
                // z: <x|z|xz> + 2/root3<z|z|z2> + <y|z|yz>
                dx = p[k + 12] + p[k + 19] - p[k + 23] * xt + p[k + 39];
                dy = -p[k + 13] - p[k + 24] * xt + p[k + 32] + p[k + 38];
                dz = p[k + 17] + 2.0 * p[k + 25] * xt + p[k + 31];
                dip[0][1] -= dx * hyfpd;
                dip[1][1] -= dy * hyfpd;
                dip[2][1] -= dz * hyfpd;
            }
        }
        for (j = 1; j <= 3; ++j)
            dip[j - 1][0] += fpc_8 * fpc_1 * 1e-10 * q[i] * coord[j-1][i];
    }
    for (j = 0; j < 3; ++j) dip[j][2] = dip[j][1] + dip[j][0];
    for (j = 0; j < 3; ++j)
        dip[3][j] = std::sqrt(dip[0][j] * dip[0][j] + dip[1][j] * dip[1][j] +
                              dip[2][j] * dip[2][j]);
    if (force) {
        dipvec[1] = dip[0][2];
        dipvec[2] = dip[1][2];
        dipvec[3] = dip[2][2];
    }
    if (chargd) {
        // RESET ION'S POSITION AGAIN, SO THAT IT IS BACK WHERE IT STARTED
        for (i = 1; i <= 3; ++i)
            for (j = 1; j <= numat; ++j) coord[i-1][j] += center[i];
    }
    if (mode == 0) return dip[3][2];
    if (iw > 0) {
        char buf[256];
        std::snprintf(buf, sizeof(buf),
                      " DIPOLE       X         Y         Z       TOTAL\n"
                      " POINT-CHG. %9.3f%9.3f%9.3f%9.3f\n"
                      " HYBRID   %9.3f%9.3f%9.3f%9.3f\n"
                      " SUM      %9.3f%9.3f%9.3f%9.3f",
                      dip[0][0], dip[1][0], dip[2][0], dip[3][0],
                      dip[0][1], dip[1][1], dip[2][1], dip[3][1],
                      dip[0][2], dip[1][2], dip[2][2], dip[3][2]);
        to_screen(buf);
    }
    return dip[3][2];
}
