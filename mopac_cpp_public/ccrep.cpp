// ccrep.cpp — C++ translation of MOPAC 2016 "ccrep.F90".

#include "ccrep.h"

#include <algorithm>
#include <cmath>

#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace parameters_C;

void ccrep(int ni, int nj, double& r, double gab, double& enuclr) {
    r = r * funcon_C::a0;
    double alpni = alp[ni];
    double alpnj = alp[nj];

    double enuc = tore[ni] * tore[nj] * gab;
    double fff;
    if (ni < 101 && nj < 101) fff = xfac[ni][nj];
    else fff = 0.0;

    if (molkst_C::method_pm7) {
        if (std::abs(fff) < 1e-5) {
            if (ni < 99 && nj < 99) {
                // Lazy-fill of the shared xfac/alpb tables: only written on
                // first use of the pair, but ccrep runs inside the dcart
                // OpenMP atom-pair loop -> guard the write.
#ifdef _OPENMP
#pragma omp critical(ccrep_xfac_alpb)
#endif
                {
                    xfac[ni][nj] = 0.5 * (xfac[ni][ni] + xfac[nj][nj]);
                    alpb[ni][nj] = 0.5 * (alpb[ni][ni] + alpb[nj][nj]);
                    fff = xfac[ni][nj];
                }
            } else {
                fff = 0.0;
            }
        }
    }

    double abond = 0.0;
    if (std::abs(fff) > 1e-5) {
        abond = alpb[ni][nj];
        if (abond < 1e-6) abond = 1.2;
        double scale;
        if (!molkst_C::method_mndod) {
            if (molkst_C::method_pm6 || molkst_C::method_pm7) {
                scale = 1.0 + 2.0 * fff * exp(-abond * (r + 0.0003 * r * r * r * r * r * r));
                int i = std::max(ni, nj);
                int j = std::min(ni, nj);
                if (j == 1) {
                    if (i == 6 || i == 7) {
                        scale = 1.0 + 2.0 * fff * exp(-abond * r * r);
                    } else if (i == 8) {
                        scale = 1.0 + 2.0 * fff * exp(-abond * r * r) -
                                par3 * exp(-par4 * r * 2);
                    }
                } else if (j == 6 && i == 6) {
                    scale = scale + par1 * exp(-par2 * r);
                } else if (j == 8 && i == 14) {
                    double ax = 1.0;
                    (void)ax;
                    scale = scale - 0.7e-3 * exp(-(r - 2.9) * (r - 2.9));
                }
            } else {
                if (molkst_C::method_am1 &&
                    ((ni == 42 && nj == 1) || (ni == 1 && nj == 42))) {
                    scale = 1.0 + r * 2.0 * fff * exp(-abond * r);
                } else {
                    scale = 1.0 + 2.0 * fff * exp(-abond * r);
                }
            }
        } else {
            if (ni == nj) {
                scale = 1.0 + 2.0 * exp(-abond * r);
            } else if (nj == 11 || nj == 12 || nj == 13) {
                scale = 1.0 + exp(-abond * r) + exp(-alp[ni] * r);
            } else {
                scale = 1.0 + exp(-abond * r) + exp(-alp[nj] * r);
            }
        }
        enuclr = enuc * scale;
    } else {
        double scale;
        if (molkst_C::method_pm6 || molkst_C::method_pm7) {
            if ((ni > 56 && ni < 72) || (nj > 56 && nj < 72)) {
                scale = 10.0 * exp(-3.0 * r);
            } else {
                scale = 10.0 * exp(-2.18 * r);
            }
        } else {
            double eni = exp(-alpni * r);
            double enj = exp(-alpnj * r);
            scale = eni + enj;
            int nt = ni + nj;
            if (nt == 8 || nt == 9) {
                if (ni == 7 || ni == 8) scale = scale + (r - 1.0) * eni;
                if (nj == 7 || nj == 8) scale = scale + (r - 1.0) * enj;
            }
        }
        enuclr = std::abs(scale * enuc) + enuc;
    }

    double scale = 0.0;
    int i;
    if (molkst_C::method_pm6 || molkst_C::method_pm7) {
        double ax = guess2[ni][1] * (r - guess3[ni][1]) * (r - guess3[ni][1]);
        if (ax < 25.0) scale += tore[ni] * tore[nj] / r * guess1[ni][1] * exp(-ax);
        ax = guess2[nj][1] * (r - guess3[nj][1]) * (r - guess3[nj][1]);
        if (ax < 25.0) scale += tore[ni] * tore[nj] / r * guess1[nj][1] * exp(-ax);
        i = (abond > 1e-4) ? 0 : 4;
    } else {
        i = 4;
        if (fff > 1e-4) i = 0;
    }
    for (int ig = 1; ig <= i; ++ig) {
        if (std::abs(guess1[ni][ig]) > 0.0) {
            double ax = guess2[ni][ig] * (r - guess3[ni][ig]) * (r - guess3[ni][ig]);
            if (ax <= 25.0) scale += tore[ni] * tore[nj] / r * guess1[ni][ig] * exp(-ax);
        }
        if (std::abs(guess1[nj][ig]) <= 0.0) continue;
        double ax = guess2[nj][ig] * (r - guess3[nj][ig]) * (r - guess3[nj][ig]);
        if (ax > 25.0) continue;
        scale += tore[ni] * tore[nj] / r * guess1[nj][ig] * exp(-ax);
    }
    enuclr += scale;

    if (molkst_C::method_pm6 || molkst_C::method_pm7) {
        double ax = r / (pow((double)ni, 0.3333) + pow((double)nj, 0.3333));
        if (ax < 3.0) {
            scale = 1e-8 / pow(ax, 12);
            enuclr += std::min(scale, 1e5);
        }
    }
}
