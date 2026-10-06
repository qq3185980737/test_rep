// dhc.cpp — C++ translation of MOPAC 2016 "dhc.F90".
// D/H matrix contributions to the derivative energy for an atom pair.
// p/pa/pb are passed as 0-based contiguous arrays (dhc.h ABI); internally
// all packed matrices use the F90 1-based convention via wrapper vectors.
#include "dhc.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include "molkst_C.h"
#include "h1elec.h"
#include "rotate.h"
#include "helect.h"
#include "fock2.h"
using namespace molkst_C;

// Saved state (F90 SAVE). thread_local: dcart may run the atom-pair loop
// under OpenMP, so each thread needs its own integral scratch arrays and
// one-time init flags.
static thread_local bool ignor1 = false, ignor2 = false, lpoint = false, cutoff = false;
static thread_local int icalcn = 0;
static thread_local double cutof = 0.0, wlim = 0.0;
static thread_local std::vector<double> w(2027, 0.0), wk(2027, 0.0);  // 1-based, 2026 usable

// point: MOPAC 2016 source tree contains no point.F90 / repd.F90 (solid-state
// only; dhc's lpoint branch requires id>0 and is never entered for molecular
// systems, id==0). Kept as a stub per accepted semantics.
static void point_stub(double, int, int, int&, double*, double*, double& enuclr) {
    enuclr = 0.0;
}

void dhc(double* p, double* pa, double* pb, double* xi, int* nat, int if_,
         int il, int jf, int jl, double& dener, int mode) {
    if (icalcn != molkst_C::numcal) {
        cutof = std::min(196.0, (2.0 / 3.0 * molkst_C::cutofp) * (2.0 / 3.0 * molkst_C::cutofp));
        lpoint = false;
        icalcn = molkst_C::numcal;
        ignor1 = false;
        ignor2 = false;
        wlim = (molkst_C::id == 0) ? 0.0 : 4.0;
    }
    // xi is 3x2 column-major: atom1 = xi[0..2], atom2 = xi[3..5]
    double rij = (xi[0] - xi[3]) * (xi[0] - xi[3]) +
                 (xi[1] - xi[4]) * (xi[1] - xi[4]) +
                 (xi[2] - xi[5]) * (xi[2] - xi[5]);
    if (mode == 1 && molkst_C::id != 0) {
        ignor1 = (rij > 225.0);
        ignor2 = (rij > (4.0 / 3.0 * molkst_C::cutofp) * (4.0 / 3.0 * molkst_C::cutofp));
        if (ignor1 && ignor2) { dener = 0.0; return; }
    }
    if (mode == 2 && molkst_C::id != 0 && ignor2 && ignor1) { dener = 0.0; return; }

    std::vector<int> nfirst_loc(3), nlast_loc(3);
    nfirst_loc[1] = 1;
    nlast_loc[1] = il - if_ + 1;
    nfirst_loc[2] = nlast_loc[1] + 1;
    nlast_loc[2] = nfirst_loc[2] + jl - jf;
    int linear = (nlast_loc[2] * (nlast_loc[2] + 1)) / 2;
    // wrap 0-based input arrays into 1-based vectors (used by fock2 and helect)
    std::vector<double> pv(linear + 1, 0.0), pav(linear + 1, 0.0), pbv(linear + 1, 0.0);
    for (int t = 1; t <= linear; ++t) {
        pv[t] = p[t - 1];
        pav[t] = pa[t - 1];
        pbv[t] = pb[t - 1];
    }
    std::vector<double> f(linear + 1, 0.0), h(linear + 1, 0.0);
    int kr = 0;
    int ia = nfirst_loc[2], ic = nlast_loc[2];
    int ja = nfirst_loc[1], jc = nlast_loc[1];
    int j = 1, i = 2;
    int nj = nat[0], ni = nat[1];
    double e1b[45] = {}, e2a[45] = {};
    double shmat[81] = {};
    double enuclr_loc = 0.0;

    if (!ignor1) {
        h1elec(nj, ni, &xi[0], &xi[3], shmat);
        if (nat[0] == 102 || nat[1] == 102) {
            int k = (ic * (ic + 1)) / 2;
            for (i = 1; i <= k; ++i) h[i] = 0.0;
        } else {
            int j1 = 0;
            for (j = ia; j <= ic; ++j) {
                int jj = j * (j - 1) / 2;
                j1 = j1 + 1;
                int i1 = 0;
                for (i = ja; i <= jc; ++i) {
                    jj = jj + 1;
                    i1 = i1 + 1;
                    h[jj] = shmat[(i1 - 1) + (j1 - 1) * 9];
                    f[jj] = shmat[(i1 - 1) + (j1 - 1) * 9];
                }
            }
        }
    }
    if (!ignor2) {
        if (molkst_C::id > 0 && mode == 1) lpoint = (rij > cutof);
        if (lpoint) {
            rij = std::sqrt(rij);
            point_stub(rij, ni, nj, kr, e2a, e1b, enuclr_loc);
        } else {
            kr = 1;
            rotate(ni, nj, &xi[3], &xi[0], &w[kr], kr, e2a, e1b, enuclr_loc);
        }
        if (molkst_C::id != 0) {
            for (i = 1; i <= kr - 1; ++i) wk[i] = w[i];
            if (mode == 1) cutoff = (w[1] < wlim);
            if (cutoff) for (i = 1; i <= kr - 1; ++i) wk[i] = 0.0;
        }
        int i2 = 0;
        for (int i1 = ja; i1 <= jc; ++i1) {
            int ii = i1 * (i1 - 1) / 2 + ja - 1;
            for (int j1 = ja; j1 <= i1; ++j1) {
                ii = ii + 1;
                i2 = i2 + 1;
                h[ii] += e1b[i2 - 1];
                f[ii] += e1b[i2 - 1];
            }
        }
        i2 = 0;
        for (int i1 = ia; i1 <= ic; ++i1) {
            int ii = i1 * (i1 - 1) / 2 + ia - 1;
            for (int j1 = ia; j1 <= i1; ++j1) {
                ii = ii + 1;
                i2 = i2 + 1;
                h[ii] += e2a[i2 - 1];
                f[ii] += e2a[i2 - 1];
            }
        }
        i = -2;
        fock2(f, pv, pav, w, w, wk, i, nfirst_loc, nlast_loc, 1);
    } else {
        enuclr_loc = 0.0;
    }
    double ee = helect(nlast_loc[2], pav.data(), h.data(), f.data());
    if (molkst_C::uhf) {
        for (i = 1; i <= linear; ++i) f[i] = h[i];
        i = -2;
        fock2(f, pv, pbv, w, w, wk, i, nfirst_loc, nlast_loc, 1);
        ee = ee + helect(nlast_loc[2], pbv.data(), h.data(), f.data());
    } else {
        ee = ee * 2.0;
    }
    dener = ee + enuclr_loc;
}
