// enpart.cpp — C++ translation of MOPAC 2016 "enpart.F90" (628 lines).
//
// Energy partitioning within the (U)MNDO scheme (S. Olivella, Barcelona 1979;
// extended to AM1/PM3 by J.J.P. Stewart).  Nothing is changed on exit except
// the optional ENPART(...) fragmentation printing path.
//
// Known deviations (recorded):
//  - rotate() is wired to the current C++ rotate.cpp, whose two-electron
//    kernels (rotatd_/nddo_to_point_/elenuc_) are closed-source placeholders
//    (mozy_integral_stubs.cpp).  e(n,2)=g and the e1b/e2a-derived parts of
//    e(n,3) therefore carry placeholder numerics until M03 integral kernels
//    are fully translated.
//  - The ENPART(...) fragmentation path prints and calls post_scf_corrections
//    with local set_a/set_b vectors; the fragment-aware semantics of the
//    PM6_DH corrections (global set_a/set_b) are not yet wired, so the
//    per-fragment H-bond/dispersion numbers are not fragment-scoped.
#include "enpart.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "chanel_C.h"
#include "funcon_C.h"
#include "elemts_C.h"
#include "rotate.h"
#include "reada.h"
#include "post_scf_corrections.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace common_arrays_C;
using namespace molkst_C;
using namespace parameters_C;
using namespace elemts_C;
using namespace funcon_C;

void enpart() {
    int npair = (numat * (numat + 1)) / 2;
    // Fortran 1-based arrays: ea(..,1..2), e(..,1..4), ex(..,1..3), w2(1..8100).
    std::vector<std::vector<double>> ea(npair + 1, std::vector<double>(3, 0.0));
    std::vector<std::vector<double>> e(npair + 1, std::vector<double>(5, 0.0));
    std::vector<std::vector<double>> ex(npair + 1, std::vector<double>(4, 0.0));
    std::vector<double> w2(8101, 0.0);

    if (!uhf) pb = pa;
    bool large = (keywrd.find(" LARGE") != std::string::npos);

    // ---- ONE-CENTER ENERGIES -------------------------------------------------
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii], ib = nlast[ii], ni = nat[ii];
        ea[ii][1] = 0.0;
        int norb = ib - ia + 1;
        if (norb >= 5) {                       // d orbitals
            double t = udd[ni];
            for (int j = ia + 4; j <= ib; ++j) ea[ii][1] += p[(j * (j + 1)) / 2] * t;
        }
        if (norb >= 3) {                       // p orbitals
            double t = upp[ni];
            for (int j = ia + 1; j <= ia + 3; ++j) ea[ii][1] += p[(j * (j + 1)) / 2] * t;
        }
        ea[ii][1] += p[(ia * (ia + 1)) / 2] * uss[ni];
    }
    if (large) std::printf("\n\n\n          TOTAL ENERGY PARTITIONING\n");
    if (large) std::printf("\n          ALL ENERGIES ARE IN ELECTRON VOLTS\n");

    // F90 advances kl to a multiple of 10 capped at numat (no other effect).
    int kl = 1;
    kl += 10;
    kl = std::min(kl, numat);
    while (numat > kl) { kl += 10; kl = std::min(kl, numat); }

    double eau = 0.0, eae = 0.0;
    for (int i = 1; i <= numat; ++i) { eau += ea[i][1]; eae += ea[i][2]; }
    double tone = eau + eae;

    // ---- TWO-CENTER RESONANCE ENERGY e(n,1) --------------------------------
    int n = 1;
    for (int ii = 2; ii <= numat; ++ii) {
        e[n][1] = 0.0;
        int ia = nfirst[ii], ib = nlast[ii];
        double oneii = (nat[ii] == 102) ? 0.0 : 1.0;
        for (int jj = 1; jj <= ii - 1; ++jj) {
            n = n + 1;
            int ja = nfirst[jj], jb = nlast[jj];
            double onejj = (nat[jj] == 102) ? 0.0 : 1.0;
            e[n][1] = 0.0;
            for (int i = ia; i <= ib; ++i) {
                int ka = (i * (i - 1)) / 2;
                for (int k = ja; k <= jb; ++k) {
                    int ik = ka + k;
                    e[n][1] += 2.0 * p[ik] * h[ik] * oneii * onejj;
                }
            }
        }
        n = n + 1;
    }

    // ---- CORE-CORE REPULSION e(n,2) AND CORE-ELECTRON ATTRACTION e(n,3) -----
    n = 0;
    int kk = 0;
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii], ib = nlast[ii];
        bool si = (ib < ia);
        int ni = nat[ii];
        int iss = (ia * (ia + 1)) / 2;
        for (int jj = 1; jj <= ii - 1; ++jj) {
            n = n + 1;
            e[n][2] = 0.0;
            e[n][3] = 0.0;
            int ja = nfirst[jj], jb = nlast[jj];
            bool sj = (jb < ja);
            int nj = nat[jj];
            int jss = (ja * (ja + 1)) / 2;
            double e1b[45], e2a[45];
            for (int q = 0; q < 45; ++q) { e1b[q] = 0.0; e2a[q] = 0.0; }
            double g = 0.0;
            int kr = 0;
            double xi[3] = { coord[0][ii], coord[1][ii], coord[2][ii] };
            double xj[3] = { coord[0][jj], coord[1][jj], coord[2][jj] };
            rotate(ni, nj, xi, xj, w2.data(), kr, e1b, e2a, g);
            if (w2[1] < -1.0e20) return;       // dummy use of w2 (F90 sentinel)
            if (ib >= ia && jb >= ja) {
                kk = kk + 1;
                e[n][2] = g;
                if (ib >= ia) e[n][3] += (pa[iss] + pb[iss]) * e1b[0];
                if (jb >= ja) e[n][3] += (pa[jss] + pb[jss]) * e2a[0];
            } else {                            // sparkle contribution
                e[n][2] = g;
                if (ib >= ia) e[n][3] = (pa[iss] + pb[iss]) * e1b[0];
                if (jb >= ja) e[n][3] = (pa[jss] + pb[jss]) * e2a[0];
            }
            // <ZZ|other terms> along the side of the W square (steps of 1).
            int kinc = ((jb - ja + 1) * (jb - ja + 2)) / 2 - 1;
            int jap1 = ja + 1;
            int i = 1;
            for (int k = jap1; k <= jb; ++k) {
                int kc = (k * (k - 1)) / 2;
                for (int l = ja; l <= k; ++l) {
                    int kll = kc + l;
                    double bb = (k == l) ? 1.0 : 2.0;
                    i = i + 1;
                    if (!si) kk = kk + 1;
                    e[n][3] += (pa[kll] + pb[kll]) * bb * e2a[i - 1];
                }
            }
            // <sz|ZZ> terms down the side (steps of kinc+1).
            int iap1 = ia + 1;
            int l = 1;
            for (int ii2 = iap1; ii2 <= ib; ++ii2) {
                int ka = (ii2 * (ii2 - 1)) / 2;
                for (int j = ia; j <= ii2; ++j) {
                    int ij = ka + j;
                    double aa = (ii2 == j) ? 1.0 : 2.0;
                    l = l + 1;
                    if (!sj) kk = kk + 1;
                    e[n][3] += (pa[ij] + pb[ij]) * aa * e1b[l - 1];
                    if (!sj) kk = kk + kinc;
                }
            }
        }
        // One-center coulomb and exchange terms for atom ii.
        double sum1 = 0.0;
        for (int i = ia; i <= ib; ++i) {
            double aa = 2.0;
            for (int j = ia; j <= i; ++j) {
                if (i == j) aa = 1.0;
                int jij = (i * (i - 1)) / 2 + j;
                double sum = 0.0;
                for (int k = ia; k <= ib; ++k) {
                    double bb = 2.0;
                    for (int l = ia; l <= k; ++l) {
                        if (k == l) bb = 1.0;
                        int jkl = (k * (k - 1)) / 2 + l;
                        int im = std::max(i, k), km = std::min(i, k);
                        int kik = (im * (im - 1)) / 2 + km;
                        im = std::max(i, l); km = std::min(i, l);
                        int kil = (im * (im - 1)) / 2 + km;
                        int jm = std::max(j, k); km = std::min(j, k);
                        int kjk = (jm * (jm - 1)) / 2 + km;
                        jm = std::max(j, l); km = std::min(j, l);
                        int kjl = (jm * (jm - 1)) / 2 + km;
                        kk = kk + 1;
                        sum += aa * bb * w2[kk] *
                            ((pa[jij] + pb[jij]) * (pa[jkl] + pb[jkl]) -
                             (pa[kik] * pa[kjl] + pa[kil] * pa[kjk] +
                              pb[kik] * pb[kjl] + pb[kil] * pb[kjk]) * 0.5);
                    }
                }
                sum1 += sum;
            }
        }
        ea[ii][2] = sum1 * 0.5;
        n = n + 1;
        e[n][2] = 0.0;
        e[n][3] = 0.0;
    }

    eau = 0.0; eae = 0.0;
    for (int i = 1; i <= numat; ++i) { eau += ea[i][1]; eae += ea[i][2]; }
    tone = eau + eae;

    // ---- TWO-CENTER COULOMB e(n,4) AND EXCHANGE ex(n,1) ---------------------
    n = 1;
    kk = 0;
    int norb1 = nlast[1] - nfirst[1] + 1;
    kk = ((norb1 * (norb1 + 1)) / 2) * ((norb1 * (norb1 + 1)) / 2);
    for (int ii = 2; ii <= numat; ++ii) {
        e[n][4] = 0.0;
        ex[n][1] = 0.0;
        int ia = nfirst[ii], ib = nlast[ii];
        for (int jj = 1; jj <= ii - 1; ++jj) {
            n = n + 1;
            e[n][4] = 0.0;
            ex[n][1] = 0.0;
            int ja = nfirst[jj], jb = nlast[jj];
            for (int i = ia; i <= ib; ++i) {
                int ka = (i * (i - 1)) / 2;
                for (int j = ia; j <= i; ++j) {
                    int kb = (j * (j - 1)) / 2;
                    int ij = ka + j;
                    double aa = (i == j) ? 1.0 : 2.0;
                    double pij = (pa[ij] + pb[ij]);
                    for (int k = ja; k <= jb; ++k) {
                        int kc = (k * (k - 1)) / 2;
                        int ik = ka + k;
                        int jk = kb + k;
                        for (int l = ja; l <= k; ++l) {
                            int il = ka + l;
                            int jl = kb + l;
                            int kll = kc + l;
                            double bb = (k == l) ? 1.0 : 2.0;
                            kk = kk + 1;
                            double g2 = w2[kk];
                            e[n][4] += aa * bb * g2 * pij * (pa[kll] + pb[kll]);
                            ex[n][1] -= 0.5 * aa * bb * g2 *
                                (pa[ik] * pa[jl] + pa[il] * pa[jk] +
                                 pb[ik] * pb[jl] + pb[il] * pb[jk]);
                        }
                    }
                }
            }
        }
        int norbs_here = ib - ia + 1;
        kk = kk + ((norbs_here * (norbs_here + 1)) / 2) * ((norbs_here * (norbs_here + 1)) / 2);
        n = n + 1;
    }
    int numat1 = npair;
    for (int q = 1; q <= 4; ++q) e[numat1][q] = 0.0;
    for (int q = 1; q <= 3; ++q) ex[numat1][q] = 0.0;

    // ex(,2) = resonance + exchange; diagonal holds one-center e-e repulsion.
    for (int q = 1; q <= numat1; ++q) ex[q][2] = e[q][1] + ex[q][1];
    for (int i = 1; i <= numat; ++i) ex[(i * (i + 1)) / 2][2] = ea[i][2];
    for (int i = 1; i <= numat; ++i) e[(i * (i + 1)) / 2][3] = ea[i][1];
    for (int q = 1; q <= numat1; ++q) ex[q][3] = e[q][4] + e[q][3] + e[q][2];

    // ---- LARGE: one/two-center tables ---------------------------------------
    if (large) {
        std::printf("\n\n          ONE-CENTER TERMS\n");
        std::printf("          E-E:  ELECTRON-ELECTRON REPULSION\n");
        std::printf("          E-N:  ELECTRON-NUCLEAR ATTRACTION\n");
        std::printf("\n   ATOM      E-E       E-N    (E-E + E-N)\n");
        for (int i = 1; i <= numat; ++i) {
            int j = (i * (i + 1)) / 2;
            std::printf("  %2s%4d %10.4f %10.4f %10.4f\n",
                        elemnt[nat[i]].c_str(), i, ex[j][2], e[j][3], ex[j][2] + e[j][3]);
        }
        std::printf("\n          TWO-CENTER TERMS\n");
        std::printf("          J:   RESONANCE ENERGY          E-E: ELECTRON-ELECTRON REPULSION\n");
        std::printf("          K:   EXCHANGE ENERGY           E-N: ELECTRON-NUCLEAR ATTRACTION\n");
        std::printf("                                    N-N: NUCLEAR-NUCLEAR REPULSION\n");
        std::printf("          C:   COULOMBIC INTERACTION = E-E + E-N + N-N\n");
        std::printf("          EE:  TOTAL OF ELECTRONIC AND NUCLEAR ENERGIES\n");
        std::printf("\n     ATOM          J        K       E-E       E-N      N-N      C        EE\n");
        std::printf("     PAIR\n");
        int ij = 0;
        for (int i = 1; i <= numat; ++i) {
            if (i < 6 || i == numat) {
                for (int j = 1; j <= i; ++j) {
                    ij = ij + 1;
                    if (i != j) {
                        std::printf("%2s%5d %2s%5d %9.4f %9.4f %9.4f %10.4f %9.4f %8.4f %9.4f\n",
                                    elemnt[nat[i]].c_str(), i, elemnt[nat[j]].c_str(), j,
                                    e[ij][1], ex[ij][1], e[ij][4], e[ij][3], e[ij][2],
                                    ex[ij][3], ex[ij][2] + ex[ij][3]);
                    } else {
                        std::printf("\n");
                    }
                }
            } else {
                for (int j = 1; j <= i; ++j) {
                    ij = ij + 1;
                    if (i != j) {
                        std::printf("%2s%5d %2s%5d %9.4f %9.4f %9.4f %10.4f %9.4f %8.4f %9.4f\n",
                                    elemnt[nat[i]].c_str(), i, elemnt[nat[j]].c_str(), j,
                                    e[ij][1], ex[ij][1], e[ij][4], e[ij][3], e[ij][2],
                                    ex[ij][3], ex[ij][2] + ex[ij][3]);
                    } else {
                        std::printf("\n   ATOM          J        K       E-E       E-N      N-N      C        EE\n");
                        std::printf("   PAIR\n");
                    }
                }
            }
        }
    } else {
        std::printf("\n          For more detail of energy partitioning, add keyword 'LARGE'\n");
    }

    // ---- ENPART(...) fragmentation (structure preserved; see deviations) ----
    if (keywrd.find("ENPART(") != std::string::npos) {
        std::size_t i = keywrd.find("ENPART(");
        std::size_t j = keywrd.find(")", i);
        keywrd[j] = ',';
        int k = 0;
        int parts[10];
        int nparts = 0;
        for (nparts = 1; nparts <= 10; ++nparts) {
            parts[nparts - 1] = (int)std::lround(reada(keywrd, (int)i + 1)) + k;
            k = parts[nparts - 1];
            i = keywrd.find(",", i);
            if (i == std::string::npos || (int)i >= (int)j) break;
        }
        if (k > numat) {
            std::printf("  Number of atoms in parts exceeds total number of atoms\n");
            std::printf("  Number of atoms in parts: %d\n", k);
            std::printf("  Number of atoms:         %d\n", numat);
        }
        keywrd[j] = ')';
        if (k <= numat) {
            std::printf("\n           In the energy partitioning, the system is to be split into %d parts\n", nparts);
            std::printf("          Part     No. of atoms     Atoms\n");
            k = 0;
            for (int i2 = 1; i2 <= nparts; ++i2) {
                std::printf("             %3d %12d %10d ' to'%4d\n", i2, parts[i2 - 1] - k, k + 1, parts[i2 - 1]);
                k = parts[i2 - 1];
            }
            std::vector<int> set_a(numat + 1, 0), set_b(numat + 1, 0);
            std::printf("                                     Contribution from:    H-bonds    +  Dispersion   =   Total\n");
            int il = 1;
            for (int nparts2 = 1; nparts2 <= nparts; ++nparts2) {
                int iu = parts[nparts2 - 1];
                itemp_1 = -99;
                std::fill(set_a.begin(), set_a.end(), 0);
                std::fill(set_b.begin(), set_b.end(), 0);
                int j2 = 0;
                for (int i2 = il; i2 <= iu; ++i2) { j2 = j2 + 1; set_a[j2] = i2; set_b[j2] = i2; }
                numcal = numcal + 1;
                double sum = 0.0;
                post_scf_corrections(sum, false);
                sum = 0.0;
                for (int i2 = il; i2 <= iu; ++i2) {
                    int k2 = (i2 * (i2 - 1)) / 2 + il - 1;
                    for (int j3 = il; j3 <= i2 - 1; ++j3) { k2 = k2 + 1; sum += ex[k2][2] + ex[k2][3]; }
                    k2 = k2 + 1;
                    sum += ex[k2][2] + e[k2][3];
                }
                std::printf(" Part %d self - energy: %14.4f eV = %13.3f Kcal/mol %10.3f %16.3f %15.3f Kcal/mol\n",
                            nparts2, sum, sum * fpc_9, E_hb, E_disp, E_disp + E_hb);
                int jl = 1;
                for (int l = 1; l <= nparts - 1; ++l) {
                    int ju = parts[l - 1];
                    sum = 0.0;
                    for (int i2 = il; i2 <= iu; ++i2) {
                        int k2 = (i2 * (i2 - 1)) / 2 + jl - 1;
                        for (int j3 = jl; j3 <= ju; ++j3) { k2 = k2 + 1; sum += ex[k2][2] + ex[k2][3]; }
                    }
                    itemp_1 = -99;
                    std::fill(set_a.begin(), set_a.end(), 0);
                    std::fill(set_b.begin(), set_b.end(), 0);
                    j2 = 0;
                    for (int i2 = il; i2 <= iu; ++i2) { j2 = j2 + 1; set_a[j2] = i2; }
                    j2 = 0;
                    for (int i2 = jl; i2 <= ju; ++i2) { j2 = j2 + 1; set_b[j2] = i2; }
                    numcal = numcal + 1;
                    double aa = 0.0;
                    post_scf_corrections(aa, false);
                    std::printf(" Parts %d - %d interaction: %11.4f eV = %13.3f Kcal/mol %10.3f %16.3f %15.3f Kcal/mol\n",
                                nparts, l, sum, sum * fpc_9, E_hb, E_disp, E_disp + E_hb);
                    jl = ju + 1;
                }
                il = iu + 1;
            }
        }
    }
    itemp_1 = 0;
    numcal = numcal + 1;
    double sum = 0.0;
    post_scf_corrections(sum, false);
    std::printf("\n                    Total contribution from hydrogen bonds:    %10.3f disp.: %9.3f Tot: %10.3f kcal/mol\n",
                E_hb, E_disp, E_disp + E_hb);

    // ---- TOTALS -----------------------------------------------------------------
    double eabr = 0.0, eabx = 0.0, eabee = 0.0, eaben = 0.0, eabnn = 0.0;
    for (int i = 1; i <= numat1; ++i) {
        eabr += e[i][1];
        eabx += ex[i][1];
        eabee += e[i][4];
        eaben += e[i][3];
        eabnn += e[i][2];
    }
    double eabrx = eabr + eabx;
    double eabe = eabee + eaben + eabnn;
    double ttwo = eabrx + eabe;
    double et = tone + ttwo;
    std::printf("\n\n\n***  SUMMARY OF ENERGY PARTITION  ***\n");
    std::printf(" ---------------------------------------\n");
    std::printf("     ONE-CENTER TERMS\n");
    std::printf("\n ELECTRON-NUCLEAR  (ONE-ELECTRON) %17.4f EV\n", eau);
    std::printf(" ELECTRON-ELECTRON (TWO-ELECTRON) %17.4f EV\n", eae);
    std::printf("\n TOTAL OF ONE-CENTER TERMS %15.4f EV\n", tone);
    std::printf(" ---------------------------------------\n");
    std::printf("     TWO-CENTER TERMS\n");
    std::printf("\n RESONANCE ENERGY %19.4f EV\n", eabr);
    std::printf(" EXCHANGE ENERGY %19.4f EV\n", eabx);
    std::printf("\n EXCHANGE + RESONANCE ENERGY: %15.4f EV\n", eabrx);
    std::printf("\n ELECTRON-ELECTRON REPULSION %23.4f EV\n", eabee);
    std::printf(" ELECTRON-NUCLEAR ATTRACTION %23.4f EV\n", eaben);
    std::printf(" NUCLEAR-NUCLEAR REPULSION  %23.4f EV\n", eabnn);
    std::printf("\n TOTAL ELECTROSTATIC INTERACTION %18.4f EV\n", eabe);
    std::printf("\n GRAND TOTAL OF TWO-CENTER TERMS %19.4f EV\n", ttwo);
    std::printf(" ---------------------------------------\n");
    std::printf(" ETOT (EONE + ETWO) %17.4f EV\n\n", et);
}
