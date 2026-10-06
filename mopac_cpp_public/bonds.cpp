// bonds.cpp — C++ translation of MOPAC 2016 "bonds.F90".
// Bond orders and valencies. vecprt / to_screen / bonds_for_MOZYME are stubs.

#include "bonds.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

namespace {
void vecprt(const std::vector<double>&, int) {}
void to_screen(const char*) {}
void bonds_for_MOZYME() {}
}

void bonds() {
    if (mozyme) {
        bonds_for_MOZYME();
        return;
    }
    int i = (numat * (numat + 1)) / 2;
    int j = (norbs * (norbs + 1)) / 2;
    bondab.assign(i + 1, 0.0);
    std::vector<std::vector<double>> b(norbs + 1, std::vector<double>(norbs + 1, 0.0));
    std::vector<double> sdm(j + 1, 0.0);
    std::vector<std::vector<double>> aux(numat + 1, std::vector<double>(numat + 1, 0.0));
    std::vector<double> pspin(j + 1, 0.0), spinab(i + 1, 0.0);
    std::vector<double> spna(numat + 1), v(numat + 1), fv(numat + 1), sq(numat + 1),
        aq(numat + 1), tq(numat + 1), pm(numat + 1), sp(numat + 1), sd(numat + 1),
        spsa(numat + 1), spsq(numat + 1);

    bool ci = (keywrd.find("C.I.") != std::string::npos ||
               keywrd.find("MECI") != std::string::npos);
    bool kci = (keywrd.find("MICROS") == std::string::npos);
    bool nci = (keywrd.find("ROOT") == std::string::npos &&
                keywrd.find("OPEN") == std::string::npos);
    int nopn = nopen - nclose;
    nelecs = nclose + nclose + nopn + nalpha + nbeta;

    // Packed lower triangle -> full symmetric b.
    int k = 0;
    for (i = 1; i <= norbs; ++i)
        for (j = 1; j <= i; ++j) {
            ++k;
            b[i][j] = p[k];
            b[j][i] = p[k];
        }

    // kappa factor for UHF / ROHF.
    double zkappa;
    if (keywrd.find("UHF") != std::string::npos) {
        zkappa = 0.0;
        for (int n = 1; n <= nalpha; ++n)
            for (int m = 1; m <= nbeta; ++m) {
                double sum = 0.0;
                for (int mu = 1; mu <= norbs; ++mu) sum += c[mu][n] * cb[mu][m];
                zkappa += sum * sum;
            }
        zkappa = 1.0 / (zkappa / (double)(nalpha + nbeta) + 0.5);
    } else {
        if (!ci && nopn == 0 && nci && kci) {
            zkappa = 1.0;
        } else {
            zkappa = 1.0 / (1.0 - ((double)nopn / (double)nelecs) / 2.0);
        }
    }

    int ij = 0;
    int l = 0;
    int ll = 0;
    for (i = 1; i <= numat; ++i) {
        double a = 0.0;
        l = nfirst[i]; ll = nlast[i];
        for (j = 1; j <= i; ++j) {
            ++ij;
            int kk = nlast[j];
            int kf = nfirst[j];
            double x = 0.0;
            for (int il = l; il <= ll; ++il)
                for (int ih = kf; ih <= kk; ++ih) x += b[il][ih] * b[il][ih];
            bondab[ij] = x;
        }
        double x = -bondab[ij];
        for (j = l; j <= ll; ++j) {
            a += b[j][j];
            x += 2.0 * b[j][j];
        }
        v[i] = x;
        sd[i] = a;
    }

    k = 0;
    for (i = 1; i <= numat; ++i) {
        for (j = 1; j <= i; ++j) {
            ++k;
            bondab[k] *= zkappa;
            aux[i][j] = bondab[k];
            aux[j][i] = bondab[k];
        }
        bondab[k] = v[i];
    }
    for (i = 1; i <= numat; ++i) {
        double da = 0.0;
        for (j = 1; j <= numat; ++j) {
            if (j == i) continue;
            da += aux[i][j];
        }
        aq[i] = da;
        sq[i] = (aux[i][i] - da) / 2.0;
        fv[i] = v[i] - aq[i];
        tq[i] = aq[i] + sq[i];
        pm[i] = sd[i] - parameters_C::tore[nat[i]];
        sp[i] = tq[i] - parameters_C::tore[nat[i]];
    }

    std::printf("\n\n\n BOND ORDERS AND VALENCIES\n");
    vecprt(bondab, numat);
    to_screen("To_file: Bonds");
    if (keywrd.find(" LARGE") == std::string::npos) return;

    // Spin population analysis.
    std::printf("\n\n SELF-Q / ACTIV-Q / TOTAL-Q table\n");
    if (keywrd.find("UHF") == std::string::npos) {
        if (!ci && nopn == 0 && nci && kci) {
            std::printf("CLOSED SHELL\n");
            return;
        }
        dopen(c, norbs, norbs, nclose, nopen, fract, sdm);
        pspin = sdm;
        std::printf("ROHF\n");
    } else {
        std::printf("UHF\n");
        pspin = pa;
        for (int idx = 1; idx < (int)pspin.size(); ++idx) pspin[idx] -= pb[idx];
    }
    double sum = 0.0;
    l = 0;
    for (i = 1; i <= norbs; ++i)
        for (j = 1; j <= i; ++j) {
            double aa = 2.0;
            if (i == j) aa = 1.0;
            ++l;
            sum += aa * (pspin[l] * p[l]);
        }
    std::printf("NALPHA-NBETA= %10.5f\n", sum);

    ij = 0;
    for (i = 1; i <= numat; ++i) {
        l = nfirst[i]; ll = nlast[i];
        for (j = 1; j <= i; ++j) {
            ++ij;
            int kf = nfirst[j], kk = nlast[j];
            double x = 0.0;
            for (int il = l; il <= ll; ++il)
                for (int ih = kf; ih <= kk; ++ih) {
                    int ilih = (il >= ih) ? (il * (il - 1)) / 2 + ih
                                          : (ih * (ih - 1)) / 2 + il;
                    x += b[il][ih] * pspin[ilih];
                }
            spinab[ij] = x;
        }
    }
    k = 0;
    for (i = 1; i <= numat; ++i)
        for (j = 1; j <= i; ++j) {
            ++k;
            aux[i][j] = spinab[k];
            aux[j][i] = spinab[k];
        }
    for (i = 1; i <= numat; ++i) {
        double da = 0.0;
        for (j = 1; j <= numat; ++j) {
            if (j == i) continue;
            da += aux[i][j];
        }
        spsa[i] = da;
        spsq[i] = aux[i][i];
        spna[i] = da + aux[i][i];
    }
    vecprt(spinab, numat);
}

void dopen(const std::vector<std::vector<double>>& eigvec, int mdim, int nor,
           int ndubl, int nsingl, double frac_in, std::vector<double>& sdm) {
    (void)mdim;
    int nl1 = ndubl + 1;
    int nu1 = nsingl;
    int l = 0;
    for (int i = 1; i <= nor; ++i)
        for (int j = 1; j <= i; ++j) {
            ++l;
            double sum1 = 0.0;
            for (int n = nl1; n <= nu1; ++n) sum1 += eigvec[i][n] * eigvec[j][n];
            sdm[l] = sum1 * frac_in;
        }
}
