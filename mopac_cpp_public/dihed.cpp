// dihed.cpp — C++ translation of MOPAC 2016 "dihed.F90".
// Dihedral angle i-j-k-l (radians). For periodic systems (id != 0) the
// minimum-image vectors across the supercell are used.
// Coordinate convention: 1-based rows (coord[1]=x, coord[2]=y, coord[3]=z).

#include "dihed.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "dang.h"
#include "funcon_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace funcon_C;
using namespace common_arrays_C;
using namespace molkst_C;

void dihed(const std::vector<std::vector<double>>& xyz,
           int i, int j, int k, int l, double& angle) {
    double xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1;
    if (id == 0) {
        xi1 = xyz[0][i] - xyz[0][k];
        xj1 = xyz[0][j] - xyz[0][k];
        xl1 = xyz[0][l] - xyz[0][k];
        yi1 = xyz[1][i] - xyz[1][k];
        yj1 = xyz[1][j] - xyz[1][k];
        yl1 = xyz[1][l] - xyz[1][k];
        zi1 = xyz[2][i] - xyz[2][k];
        zj1 = xyz[2][j] - xyz[2][k];
        zl1 = xyz[2][l] - xyz[2][k];
    } else {
        // Minimum-image vectors across lattice images.
        double rmin_ik = 1.0e8, rmin_jk = 1.0e8, rmin_lk = 1.0e8;
        double Vab_ik[3] = {0.0, 0.0, 0.0}, Vab_jk[3] = {0.0, 0.0, 0.0},
               Vab_lk[3] = {0.0, 0.0, 0.0};
        for (int ii = -l11; ii <= l11; ++ii) {
            for (int jj = -l21; jj <= l21; ++jj) {
                for (int kk = -l31; kk <= l31; ++kk) {
                    double Vab[3];
                    Vab[0] = xyz[0][i] - xyz[0][k] + tvec[1][1] * ii +
                             tvec[1][2] * jj + tvec[1][3] * kk;
                    Vab[1] = xyz[1][i] - xyz[1][k] + tvec[2][1] * ii +
                             tvec[2][2] * jj + tvec[2][3] * kk;
                    Vab[2] = xyz[2][i] - xyz[2][k] + tvec[3][1] * ii +
                             tvec[3][2] * jj + tvec[3][3] * kk;
                    double Rab = Vab[0] * Vab[0] + Vab[1] * Vab[1] +
                                 Vab[2] * Vab[2];
                    if (Rab < rmin_ik) {
                        rmin_ik = Rab;
                        Vab_ik[0] = Vab[0]; Vab_ik[1] = Vab[1];
                        Vab_ik[2] = Vab[2];
                    }
                    Vab[0] = xyz[0][j] - xyz[0][k] + tvec[1][1] * ii +
                             tvec[1][2] * jj + tvec[1][3] * kk;
                    Vab[1] = xyz[1][j] - xyz[1][k] + tvec[2][1] * ii +
                             tvec[2][2] * jj + tvec[2][3] * kk;
                    Vab[2] = xyz[2][j] - xyz[2][k] + tvec[3][1] * ii +
                             tvec[3][2] * jj + tvec[3][3] * kk;
                    Rab = Vab[0] * Vab[0] + Vab[1] * Vab[1] + Vab[2] * Vab[2];
                    if (Rab < rmin_jk) {
                        rmin_jk = Rab;
                        Vab_jk[0] = Vab[0]; Vab_jk[1] = Vab[1];
                        Vab_jk[2] = Vab[2];
                    }
                    Vab[0] = xyz[0][l] - xyz[0][k] + tvec[1][1] * ii +
                             tvec[1][2] * jj + tvec[1][3] * kk;
                    Vab[1] = xyz[1][l] - xyz[1][k] + tvec[2][1] * ii +
                             tvec[2][2] * jj + tvec[2][3] * kk;
                    Vab[2] = xyz[2][l] - xyz[2][k] + tvec[3][1] * ii +
                             tvec[3][2] * jj + tvec[3][3] * kk;
                    Rab = Vab[0] * Vab[0] + Vab[1] * Vab[1] + Vab[2] * Vab[2];
                    if (Rab < rmin_lk) {
                        rmin_lk = Rab;
                        Vab_lk[0] = Vab[0]; Vab_lk[1] = Vab[1];
                        Vab_lk[2] = Vab[2];
                    }
                }
            }
        }
        xi1 = Vab_ik[0]; xj1 = Vab_jk[0]; xl1 = Vab_lk[0];
        yi1 = Vab_ik[1]; yj1 = Vab_jk[1]; yl1 = Vab_lk[1];
        zi1 = Vab_ik[2]; zj1 = Vab_jk[2]; zl1 = Vab_lk[2];
    }
    // Rotate around Z axis to put KJ along Y axis.
    double dist = std::sqrt(xj1 * xj1 + yj1 * yj1 + zj1 * zj1);
    double cosa = zj1 / dist;
    cosa = std::min(1.0, cosa);
    cosa = std::max(-1.0, cosa);
    double ddd = 1.0 - cosa * cosa;
    double xi2, xl2, yi2, yl2, costh, sinth;
    double yxdist = dist * std::sqrt(ddd);
    if (ddd <= 0.0 || yxdist <= 1.0e-6) {
        xi2 = xi1;
        xl2 = xl1;
        yi2 = yi1;
        yl2 = yl1;
        costh = cosa;
        sinth = 0.0;
    } else {
        double cosph = yj1 / yxdist;
        double sinph = xj1 / yxdist;
        xi2 = xi1 * cosph - yi1 * sinph;
        xl2 = xl1 * cosph - yl1 * sinph;
        yi2 = xi1 * sinph + yi1 * cosph;
        double yj2 = xj1 * sinph + yj1 * cosph;
        yl2 = xl1 * sinph + yl1 * cosph;
        // Rotate KJ around X axis so KJ lies along Z axis.
        costh = cosa;
        sinth = yj2 / dist;
    }
    double yi3 = yi2 * costh - zi1 * sinth;
    double yl3 = yl2 * costh - zl1 * sinth;
    dang(xl2, yl3, xi2, yi3, angle);
    if (angle < 0.0) angle = pi * 2.0 + angle;
    if (angle >= 6.28318530717959) angle = 0.0;
}
