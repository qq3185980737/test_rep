// outer1.cpp
#include "outer1.h"
#include <cmath>
#include <algorithm>
#include "molkst_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace parameters_C;
using namespace funcon_C;
using namespace common_arrays_C;

// trunk(r): smooth truncation of interatomic distance for Madelung sums
// (MOPAC 2016 solrot.F90, lines 183-239).  save semantics -> file-static.
static int icalcn_trk = 0;
static double bound1_trk = 0, bound2_trk = 0, c_trk = 0, clim_trk = 0,
              cr_trk = 0, cr2_trk = 0, range_trk = 0;
static double trunk_cpp(double r) {
    if (icalcn_trk != numcal) {
        clower = std::min(cutofp * 2.0 / 3.0, 13.0);
        bound1_trk = clower / cutofp;
        cupper = cutofp;
        bound2_trk = cupper / cutofp;
        range_trk = bound2_trk - bound1_trk;
        c_trk = -0.5 * bound1_trk * bound1_trk * cutofp / range_trk;
        cr_trk = 1.0 + bound1_trk / range_trk;
        cr2_trk = -1.0 / (cutofp * 2.0 * range_trk);
        clim_trk = c_trk + cr_trk * cupper + cr2_trk * cupper * cupper;
        icalcn_trk = numcal;
    }
    if (r > clower) {
        if (r > cupper) return clim_trk;
        return c_trk + cr_trk * r + cr2_trk * r * r;
    }
    return r;
}

// to_point(r, point, const): NDDO->point-charge smooth transition
// (MOPAC 2016 mndod.F90, lines 3501-3525).
void to_point(double rij, double& point, double& cnst) {
    double& const_ = cnst;
    point = ev * a0 / rij;
    if (rij < trunc_1)
        const_ = 1.0 - std::exp(-(rij - trunc_1) * (rij - trunc_1) * trunc_2);
    else
        const_ = 0.0;
}
void outer1(int ni_loc, int nj, const double* c1, const double* c2,
             double* w, int& kr, double* e1b, double* e2a,
             double& enuc, int mode, bool direct) {
    if (mode == 0) {
        for (int i = 0; i < 45; ++i) { e1b[i] = 0.0; e2a[i] = 0.0; }
        double rij = std::sqrt((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));
        double r = rij / a0;
        double aee = 0.5/am[ni_loc] + 0.5/am[nj];
        aee = aee * aee;
        double ww = ev / std::sqrt(r*r + aee);
        if (method_pm7) {
            double point, const_;
            to_point(rij, point, const_);
            ww = ww*const_ + (1.0-const_)*point;
        }
        int idx[9] = {0,2,5,9,14,20,27,35,44};
        for (int k = 0; k < 9; ++k) {
            e1b[idx[k]] = -ww * tore[nj];
            e2a[idx[k]] = -ww * tore[ni_loc];
        }
        enuc = tore[ni_loc] * tore[nj] * ww;
        if (!direct) {
            w[0] = ww;
            if (natorb[ni_loc]*natorb[nj] != 0) kr++;
        }
    } else {
        w[0] = 0.0;
        // F90 clears e1b/e2a(1:10) only in solid-state mode
        for (int i = 0; i < 10; ++i) { e1b[i] = 0.0; e2a[i] = 0.0; }
        for (int i = -l1u; i <= l1u; ++i)
        for (int j = -l2u; j <= l2u; ++j)
        for (int k = -l3u; k <= l3u; ++k) {
            double c2s[3];
            for (int l = 0; l < 3; ++l)
                c2s[l] = c2[l] + tvec[l][1]*i + tvec[l][2]*j + tvec[l][3]*k - c1[l];
            double r = std::sqrt(c2s[0]*c2s[0] + c2s[1]*c2s[1] + c2s[2]*c2s[2]);
            r = trunk_cpp(r);
            w[0] += a0 * ev / r;
        }
        int idx[9] = {0,2,5,9,14,20,25,35,44};
        for (int k = 0; k < 9; ++k) {
            e1b[idx[k]] = -w[0] * tore[nj];
            e2a[idx[k]] = -w[0] * tore[ni_loc];
        }
        enuc = tore[ni_loc] * tore[nj] * w[0];
        if (natorb[ni_loc]*natorb[nj] != 0) kr++;
    }
}
