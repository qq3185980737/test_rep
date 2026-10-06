// anavib.cpp — C++ translation of MOPAC 2016 "anavib.F90".
#pragma warning(disable: 4459)
// Description of vibrations: interatomic distances, energy contribution of
// each atom pair, and formatted report. Output unit iw maps to stdout.

#include "anavib.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "to_screen_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace funcon_C;

namespace {
// vanrad(107): hard-coded van der Waals radii (Angstroms).
const std::vector<double> vanrad = {
    0.32, 0.93, 1.23, 0.90, 0.82, 0.77, 0.75, 0.73, 0.72, 0.71,
    1.54, 1.36, 1.18, 1.11, 1.06, 1.02, 0.99, 0.98,
    2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25,
    1.26, 1.22, 1.20, 1.16, 1.14, 1.12,
    2.16, 1.91, 1.62, 1.45, 1.34, 1.30, 1.27, 1.25, 1.25, 1.28, 1.34, 1.48, 1.44, 1.41,
    1.40, 1.36, 1.33, 1.31,
    2.35, 1.98, 1.69, 1.65, 1.65, 1.64, 1.63, 1.62, 1.85, 1.61, 1.59, 1.59, 1.58, 1.57, 1.56, 1.56, 1.56,
    1.44, 1.34, 1.30, 1.28, 1.26, 1.27, 1.30, 1.34, 1.49, 1.48, 1.47, 1.46, 1.46, 1.45, 1.45,
    1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45,
    1.45, 1.45, 1.45, 1.45, 1.45};
}

void anavib(std::vector<double>& eigs, std::vector<double>& dipt, int n3,
            std::vector<std::vector<double>>& vibs, std::vector<double>& rij,
            int nv, std::vector<double>& hess, std::vector<double>& f) {
    std::vector<int> ijf(11, 0);
    std::vector<double> fij(11, 0.0);
    (void)n3;
    int l, i, j, iline, k, j3, linear, j1, i1, jj = 0, ii, ij, j2, i3, i2;

    // COMPUTE INTERATOMIC DISTANCES.
    l = 0;
    for (i = 1; i <= molkst_C::numat; ++i) {
        for (j = 1; j <= i - 1; ++j) {
            ++l;
            rij[l] = std::sqrt((coord[0][j] - coord[0][i]) * (coord[0][j] - coord[0][i]) +
                               (coord[1][j] - coord[1][i]) * (coord[1][j] - coord[1][i]) +
                               (coord[2][j] - coord[2][i]) * (coord[2][j] - coord[2][i])) +
                     1.0e-10;
        }
    }

    // ANALYSE VIBRATIONS
    std::printf("\n\n          DESCRIPTION OF VIBRATIONS\n\n");
    iline = 0;
    for (k = 1; k <= nv; ++k) {
        bool vib1 = true, vib2 = true, vib3 = true, vib4 = true, vib5 = true;
        j3 = 0;
        l = 0;
        double tot = 0.0;
        linear = 0;
        j1 = -2;
        for (j = 1; j <= molkst_C::numat; ++j) {
            j1 = j1 + 3;
            i1 = -2;
            for (i = 1; i <= j - 1; ++i) {
                i1 = i1 + 3;
                double vdw = (vanrad[nat[i]] + vanrad[nat[j]]) * 1.5;
                ++l;
                f[l] = 0.0;
                if (rij[l] >= vdw) continue;

                // CALCULATE ENERGY TERM BETWEEN THE TWO ATOMS
                double eab = 0.0;
                for (jj = j1; jj <= j1 + 2; ++jj) {
                    for (ii = i1; ii <= i1 + 2; ++ii) {
                        eab += vibs[jj][k] * hess[(jj * (jj - 1)) / 2 + ii] * vibs[ii][k];
                    }
                }
                double eb = 0.0;
                for (jj = j1; jj <= j1 + 2; ++jj) {
                    for (ii = j1; ii <= jj; ++ii) {
                        eb += vibs[jj][k] * hess[(jj * (jj - 1)) / 2 + ii] * vibs[ii][k] * 2.0;
                    }
                    eb -= vibs[jj][k] * hess[(jj * (jj + 1)) / 2] * vibs[jj][k];
                }
                double ea = 0.0;
                for (jj = i1; jj <= i1 + 2; ++jj) {
                    for (ii = i1; ii <= jj; ++ii) {
                        ea += vibs[jj][k] * hess[(jj * (jj - 1)) / 2 + ii] * vibs[ii][k] * 2.0;
                    }
                    ea -= vibs[jj][k] * hess[(jj * (jj + 1)) / 2] * vibs[jj][k];
                }
                ++linear;
                f[l] = ea + eab * 2.0 + eb;
                tot += f[l];
            }
        }
        if (k == nv || symmetry_C::jndex[k] != symmetry_C::jndex[k + 1] ||
            symmetry_C::namo[k] != symmetry_C::namo[k + 1]) {
            if (std::fabs(tot) >= 1.0e-5) {
                // NOW TO SORT F INTO DECREASING ORDER
                for (i = 1; i <= 10; ++i) {
                    double sum = -100.0;
                    for (j = 1; j <= l; ++j) {
                        if (std::fabs(f[j]) <= sum) continue;
                        jj = j;
                        sum = std::fabs(f[j]);
                    }
                    if (sum < 0.0) goto sort_done;
                    fij[i] = sum;
                    f[jj] = -1.0e-9;
                    ijf[i] = jj;
                }
                i = 10;
sort_done:
                linear = i;
                double scl = 1.0 / tot;
                for (ij = 1; ij <= linear; ++ij) {
                    j = (int)(0.5 * (0.99 + std::sqrt(1.0 + 8.0 * ijf[ij])));
                    i = ijf[ij] - (j * (j - 1)) / 2;
                    j = j + 1;
                    double xj = coord[0][j], yj = coord[1][j], zj = coord[2][j];
                    j1 = 3 * j - 2;
                    j2 = j1 + 1;
                    j3 = j2 + 1;
                    i3 = 0;
                    double xi = coord[0][i], yi = coord[1][i], zi = coord[2][i];
                    i1 = 3 * i - 2;
                    i2 = i1 + 1;
                    i3 = i2 + 1;
                    double x = vibs[j1][k] - vibs[i1][k];
                    double y = vibs[j2][k] - vibs[i2][k];
                    double z = vibs[j3][k] - vibs[i3][k];
                    double e = fij[ij] * scl * 100.0;
                    double shift = x * x + y * y + z * z + 1.0e-30;
                    if (!(std::fabs(e) > 10.0 || (ij < 5 && std::fabs(e) > 0.1))) continue;
                    shift = std::sqrt(shift);
                    double radial =
                        std::pow((x * (xi - xj) + y * (yi - yj) + z * (zi - zj)) /
                                     (shift * rij[ijf[ij]]),
                                 2.0) *
                        100.0;
                    double ans = std::min(999.0, std::max(-99.9,
                        100.0 * std::sqrt(fij[ij] * 1.0e5 * fpc_10) / (fpc_8 * pi * 2.0) /
                            eigs[k]));
                    std::string si = elemnt[nat[i]], sj = elemnt[nat[j]];
                    if (vib1) {
                        std::printf("\n VIBRATION%5d%5d%4s    ATOM PAIR        ENERGY CONTRIBUTION    RADIAL\n",
                                    k, symmetry_C::jndex[k], symmetry_C::namo[k].c_str());
                        std::printf(" FREQ.   %10.2f       %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    eigs[k], si.c_str(), i, sj.c_str(), j, e, ans, radial);
                        vib1 = false;
                    } else if (vib2) {
                        vib2 = false;
                        std::printf(" T-DIPOLE%10.4f       %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    dipt[k], si.c_str(), i, sj.c_str(), j, e, ans, radial);
                    } else if (vib3) {
                        vib3 = false;
                        std::printf(" TRAVEL  %10.4f       %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    to_screen_C::travel[k], si.c_str(), i, sj.c_str(), j, e, ans, radial);
                    } else if (vib4) {
                        vib4 = false;
                        std::printf(" RED. MASS%9.4f       %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    to_screen_C::redmas[k][1], si.c_str(), i, sj.c_str(), j, e, ans, radial);
                    } else if (vib5) {
                        vib5 = false;
                        std::printf(" EFF. MASS%9.4f       %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    to_screen_C::redmas[k][2], si.c_str(), i, sj.c_str(), j, e, ans, radial);
                    } else {
                        ++iline;
                        std::printf("                          %s%3d -- %s%3d       %9.1f%% (%5.1f%%)%8.1f%%\n",
                                    si.c_str(), i, sj.c_str(), j, e, ans, radial);
                    }
                }
            }
            if (vib1) {
                std::printf("\n VIBRATION%4d\n", k);
                std::printf(" FREQ.    %8.2f\n", eigs[k]);
            }
            if (vib2) std::printf(" T-DIPOLE %9.4f\n", dipt[k]);
            if (vib3) std::printf(" TRAVEL   %9.4f\n", to_screen_C::travel[k]);
            if (vib4) std::printf(" RED. MASS%9.4f\n", to_screen_C::redmas[k][1]);
            if (vib5) std::printf(" EFF. MASS%9.4f\n",
                std::min(9999.9999, std::max(-999.9999, to_screen_C::redmas[k][2])));
        }
    }
}
