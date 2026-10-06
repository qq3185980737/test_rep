// diat.cpp — C++ translation of MOPAC 2016 "diat.F90".
//   diat   (lines 1-154): diatomic overlap integrals between atoms ni/nj,
//   diat2  (lines 155-347): sp-only fast path (Slater overlaps),
//   ss     (lines 348-470): general STO overlap (via bfn / aff / bi tables).
// 1-based Fortran indexing is kept for di (di[0] row/col is padding).
#include "diat.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "bfn.h"
#include "coe.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "set.h"

using namespace overlaps_C;  // fact, cutof1, sa, sb, a, b, isp, ips

namespace {

// Fortran data ival/ 1,0,9, 1,3,8, 1,4,7, 1,2,6, 0,0,5/ (column-major (3,5)).
// ival(i,k): AO index for atom with L=i (1=s,2=p,3=d) and M_L-slot k (1..5).
const int ival[4][6] = {
    {0, 0, 0, 0, 0, 0},
    {0, 1, 1, 1, 1, 0},   // i=1 (s):    1,1,1,1,0
    {0, 0, 3, 4, 2, 0},   // i=2 (p):    0,3,4,2,0
    {0, 9, 8, 7, 6, 5}    // i=3 (d):    9,8,7,6,5
};

// diat2: inmb(17) — bond-type code per element.
const int inmb[18] = {0, 1, 0, 2, 2, 3, 4, 5, 6, 7, 0, 8, 8, 8, 9, 10, 11, 12};
// diat2: iii(78) — overlap case (1..6) for each bond-type pair.
const int iii[79] = {
    0,
    1, 2, 4, 2, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4,
    2, 4, 4, 4, 4, 4, 4, 3, 5, 5, 5, 5, 5, 5, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6,
    3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 3, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 3, 5,
    5, 5, 5, 5, 5, 6, 6, 6, 6, 6};

// F90 s(3,3,3) access (1-based i,j,k).
inline double& S(double s[3][3][3], int i, int j, int k) { return s[i - 1][j - 1][k - 1]; }
// F90 c(3,5,5) flat[75] access: c(i,k,l) -> flat[(l-1)*15+(k-1)*3+(i-1)].
inline double cval(const double* c, int i, int k, int l) {
    return c[(l - 1) * 15 + (k - 1) * 3 + (i - 1)];
}

// diat2: sp-only Slater-type overlap (F90 lines 155-347).
void diat2(int na, double esa, double epa, double r12, int nb, double esb,
           double epb, double s[3][3][3], double a0) {
    int jmax = std::max(inmb[na], inmb[nb]);
    int jmin = std::min(inmb[na], inmb[nb]);
    int nbond = (jmax * (jmax - 1)) / 2 + jmin;
    int ii = iii[nbond];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) s[i][j][k] = 0.0;
    double rab = r12 / a0;
    switch (ii) {
    case 1: {  // first row - first row
        set_fn(esa, esb, na, nb, rab, ii);
        S(s, 1, 1, 1) = 0.25 * std::sqrt(std::pow(sa * sb * rab * rab, 3)) *
                        (a[3] * b[1] - b[3] * a[1]);
        return;
    }
    case 2: {  // first row - second row
        set_fn(esa, esb, na, nb, rab, ii);
        double rab4 = std::pow(rab, 4) * 0.125;
        double w = std::sqrt(sa * sa * sa * std::pow(sb, 5)) * rab4;
        S(s, 1, 1, 1) = std::sqrt(1.0 / 3.0);
        S(s, 1, 1, 1) = w * S(s, 1, 1, 1) *
                        (a[4] * b[1] - b[4] * a[1] + a[3] * b[2] - b[3] * a[2]);
        if (na > 1) set_fn(epa, esb, na, nb, rab, ii);
        if (nb > 1) set_fn(esa, epb, na, nb, rab, ii);
        w = std::sqrt(sa * sa * sa * std::pow(sb, 5)) * rab4;
        S(s, isp, ips, 1) =
            w * (a[3] * b[1] - b[3] * a[1] + a[4] * b[2] - b[4] * a[2]);
        return;
    }
    case 3: {  // first row - third row
        set_fn(esa, esb, na, nb, rab, ii);
        double rab4 = std::pow(rab, 5) * 0.0625;
        double w = std::sqrt(sa * sa * sa * std::pow(sb, 7) / 22.5) * rab4;
        S(s, 1, 1, 1) =
            w * (a[5] * b[1] - b[5] * a[1] + 2.0 * (a[4] * b[2] - b[4] * a[2]));
        if (na > 1) set_fn(epa, esb, na, nb, rab, ii);
        if (nb > 1) set_fn(esa, epb, na, nb, rab, ii);
        w = std::sqrt(sa * sa * sa * std::pow(sb, 7) / 7.5) * rab4;
        S(s, isp, ips, 1) = w * (a[4] * (b[1] + b[3]) - b[4] * (a[1] + a[3]) +
                                 b[2] * (a[3] + a[5]) - a[2] * (b[3] + b[5]));
        return;
    }
    case 4: {  // second row - second row
        set_fn(esa, esb, na, nb, rab, ii);
        double rab4 = std::pow(rab, 5) * 0.0625;
        double w = std::sqrt(std::pow(sa * sb, 5)) * rab4;
        S(s, 1, 1, 1) = w * (a[5] * b[1] + b[5] * a[1] - 2.0 * a[3] * b[3]) / 3.0;
        set_fn(esa, epb, na, nb, rab, ii);
        if (na > nb) set_fn(epa, esb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 5)) * rab4;
        const double rt3 = 1.0 / std::sqrt(3.0);
        double d = a[4] * (b[1] - b[3]) - a[2] * (b[3] - b[5]);
        double e = b[4] * (a[1] - a[3]) - b[2] * (a[3] - a[5]);
        S(s, isp, ips, 1) = w * rt3 * (d + e);
        set_fn(epa, esb, na, nb, rab, ii);
        if (na > nb) set_fn(esa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 5)) * rab4;
        d = a[4] * (b[1] - b[3]) - a[2] * (b[3] - b[5]);
        e = b[4] * (a[1] - a[3]) - b[2] * (a[3] - a[5]);
        S(s, ips, isp, 1) = w * rt3 * (d - e);
        set_fn(epa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 5)) * rab4;
        S(s, 2, 2, 1) = -w * (b[3] * (a[5] + a[1]) - a[3] * (b[5] + b[1]));
        S(s, 2, 2, 2) = 0.5 * w *
                        (a[5] * (b[1] - b[3]) - b[5] * (a[1] - a[3]) -
                         a[3] * b[1] + b[3] * a[1]);
        return;
    }
    case 5: {  // second row - third row
        set_fn(esa, esb, na, nb, rab, ii);
        double rab6 = std::pow(rab, 6) * 0.03125 / std::sqrt(7.5);
        double w = std::sqrt(std::pow(sa, 5) * std::pow(sb, 7)) * rab6;
        const double rt3 = 1.0 / std::sqrt(3.0);
        S(s, 1, 1, 1) = w * (a[6] * b[1] + a[5] * b[2] -
                             2.0 * (a[4] * b[3] + a[3] * b[4]) + a[2] * b[5] +
                             a[1] * b[6]) / 3.0;
        set_fn(esa, epb, na, nb, rab, ii);
        if (na > nb) set_fn(epa, esb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa, 5) * std::pow(sb, 7)) * rab6;
        S(s, isp, ips, 1) = w * rt3 *
                            (a[6] * b[2] + a[5] * b[1] -
                             2.0 * (a[4] * b[4] + a[3] * b[3]) + a[2] * b[6] +
                             a[1] * b[5]);
        set_fn(epa, esb, na, nb, rab, ii);
        if (na > nb) set_fn(esa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa, 5) * std::pow(sb, 7)) * rab6;
        S(s, ips, isp, 1) = -w * rt3 *
                            (a[5] * (2.0 * b[3] - b[1]) - b[5] * (2.0 * a[3] - a[1]) -
                             a[2] * (b[6] - 2.0 * b[4]) + b[2] * (a[6] - 2.0 * a[4]));
        set_fn(epa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa, 5) * std::pow(sb, 7)) * rab6;
        S(s, 2, 2, 1) = -w * (b[4] * (a[1] + a[5]) - a[4] * (b[1] + b[5]) +
                              b[3] * (a[2] + a[6]) - a[3] * (b[2] + b[6]));
        S(s, 2, 2, 2) = 0.5 * w *
                        (a[6] * (b[1] - b[3]) - b[6] * (a[1] - a[3]) +
                         a[5] * (b[2] - b[4]) - b[5] * (a[2] - a[4]) -
                         a[4] * b[1] + b[4] * a[1] - a[3] * b[2] + b[3] * a[2]);
        return;
    }
    case 6: {  // third row - third row
        set_fn(esa, esb, na, nb, rab, ii);
        double rab4 = std::pow(rab, 7) / 480.0;
        double w = std::sqrt(std::pow(sa * sb, 7)) * rab4;
        const double rt3 = 1.0 / std::sqrt(3.0);
        S(s, 1, 1, 1) =
            w * (a[7] * b[1] - 3.0 * (a[5] * b[3] - a[3] * b[5]) - a[1] * b[7]) / 3.0;
        set_fn(esa, epb, na, nb, rab, ii);
        if (na > nb) set_fn(epa, esb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 7)) * rab4;
        double d = a[6] * (b[1] - b[3]) - 2.0 * a[4] * (b[3] - b[5]) +
                   a[2] * (b[5] - b[7]);
        double e = b[6] * (a[1] - a[3]) - 2.0 * b[4] * (a[3] - a[5]) +
                   b[2] * (a[5] - a[7]);
        S(s, isp, ips, 1) = w * rt3 * (d - e);
        set_fn(epa, esb, na, nb, rab, ii);
        if (na > nb) set_fn(esa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 7)) * rab4;
        d = a[6] * (b[1] - b[3]) - 2.0 * a[4] * (b[3] - b[5]) +
            a[2] * (b[5] - b[7]);
        e = b[6] * (a[1] - a[3]) - 2.0 * b[4] * (a[3] - a[5]) +
            b[2] * (a[5] - a[7]);
        S(s, ips, isp, 1) = -w * rt3 * ((-d) - e);
        set_fn(epa, epb, na, nb, rab, ii);
        w = std::sqrt(std::pow(sa * sb, 7)) * rab4;
        d = a[3] * (b[7] + b[3] + b[3]) - a[5] * (b[1] + b[5] + b[5]) -
            b[5] * a[1] + a[7] * b[3];
        S(s, 2, 2, 1) = -w * d;
        d = a[7] * (b[1] - b[3]) + b[7] * (a[1] - a[3]);
        e = a[5] * (b[5] - b[3] - b[1]) + b[5] * (a[5] - a[3] - a[1]) +
            2.0 * a[3] * b[3];
        S(s, 2, 2, 2) = 0.5 * w * (d + e);
        return;
    }
    default:
        return;  // unreachable (iii values are 1..6)
    }
}

}  // namespace

double ss_overlap(int na, int nb, int la1, int lb1, int m1,
                  double ua, double ub, double r1, double a0) {
    static thread_local std::vector<std::vector<std::vector<double>>> aff(3,
        std::vector<std::vector<double>>(3, std::vector<double>(3, 0.0)));
    static thread_local std::vector<std::vector<double>> bi(13, std::vector<double>(13, 0.0));
    static thread_local bool first = true;
    if (first) {
        first = false;
        for (int i = 0; i <= 12; ++i) { bi[i][0] = 1.0; bi[i][i] = 1.0; }
        for (int i = 0; i <= 11; ++i)
            for (int j = 1; j <= i; ++j) bi[i + 1][j] = bi[i][j] + bi[i][j - 1];
        aff[0][0][0] = 1.0;
        aff[1][0][0] = 1.0;
        aff[1][1][0] = std::sqrt(0.5);
        aff[2][0][0] = 1.5;
        aff[2][1][0] = std::sqrt(1.5);
        aff[2][2][0] = std::sqrt(0.375);
        aff[2][0][2] = -0.5;
    }
    int m = m1 - 1, lb = lb1 - 1, la = la1 - 1;
    double r = r1 / a0;
    double p = (ua + ub) * r * 0.5, b = (ua - ub) * r * 0.5;
    double quo = 1.0 / p;
    std::vector<double> af(20, 0.0), bf(20, 0.0);
    af[0] = quo * std::exp(-p);
    for (int n = 1; n <= 19; ++n) af[n] = n * quo * af[n - 1] + af[0];
    bfn(b, bf.data());  // bfn writes bf[1..13]
    double sum = 0.0;
    int lam1 = la - m, lbm1 = lb - m;
    for (int i = 0; i <= lam1; i += 2) {
        int ia = na + i - la, ic = la - i - m;
        for (int j = 0; j <= lbm1; j += 2) {
            int ib = nb + j - lb, id = lb - j - m;
            double sum1 = 0.0;
            int iab = ia + ib;
            for (int k1 = 0; k1 <= ia; ++k1)
                for (int k2 = 0; k2 <= ib; ++k2)
                    for (int k3 = 0; k3 <= ic; ++k3)
                        for (int k4 = 0; k4 <= id; ++k4)
                            for (int k5 = 0; k5 <= m; ++k5) {
                                int iaf = iab - k1 - k2 + k3 + k4 + 2 * k5;
                                for (int k6 = 0; k6 <= m; ++k6) {
                                    int ibf = k1 + k2 + k3 + k4 + 2 * k6;
                                    int sign = (1 - 2 * ((m + k2 + k4 + k5 + k6) % 2));
                                    // Fortran ss declares bf(0:19) and passes it to
                                    // bfn (which declares bf(13)): the callee's bf(1)
                                    // aliases the caller's bf(0).  Hence the C++ read
                                    // must be bf[ibf+1] (bfn wrote bf[1..13]).
                                    sum1 += bi[id][k4] * bi[ic][k3] * bi[ib][k2] *
                                            bi[ia][k1] * bi[m][k5] * bi[m][k6] * sign *
                                            af[iaf] * bf[ibf + 1];
                                }
                            }
            sum += sum1 * aff[la][m][i] * aff[lb][m][j];
        }
    }
    return sum * std::pow(r, na + nb + 1) * std::pow(ua, na) * std::pow(ub, nb) / 2.0 *
           std::sqrt(ua * ub / (fact[na + na] * fact[nb + nb]) * ((la + la + 1) * (lb + lb + 1)));
}

void diat(int ni, int nj, const std::vector<double>& xj,
          std::vector<std::vector<double>>& di) {
    static thread_local std::vector<bool> use_diat2(107, false);
    static thread_local int icalcn = -1;
    using namespace molkst_C;      // numcal
    using namespace funcon_C;      // a0
    using namespace parameters_C;  // natorb, zs, zp, zd, npq
    if (icalcn != numcal) {
        icalcn = numcal;
        std::fill(use_diat2.begin(), use_diat2.end(), false);
        for (int i = 1; i <= 17; ++i) use_diat2[i] = (natorb[i] < 5);
        use_diat2[2] = false;
        use_diat2[10] = false;
    }
    double x2 = xj[0], y2 = xj[1], z2 = xj[2];
    int pq1 = npq[ni][1], pq2 = npq[nj][1];
    for (std::size_t i = 0; i < di.size(); ++i)
        std::fill(di[i].begin(), di[i].end(), 0.0);
    double r = x2 * x2 + y2 * y2 + z2 * z2;
    if (pq1 == 0 || pq2 == 0 || r > cutof1) return;
    if (natorb[ni] == 0 || natorb[nj] == 0) return;
    std::vector<double> c(75, 0.0);
    coe(x2, y2, z2, natorb[ni], natorb[nj], c.data(), r);
    if (r < 0.001) return;
    int ia = std::min(pq1 + 1, 3);
    int ib = std::min(pq2 + 1, 3);
    int a = ia - 1, b = ib - 1;
    double s[3][3][3] = {{{0.0}}};
    if (use_diat2[ni] && use_diat2[nj]) {
        // Only s/p orbitals: run the closed diat2 path.
        diat2(ni, zs[ni], zp[ni], r, nj, zs[nj], zp[nj], s, a0);
    } else {
        double ul1[4] = {0.0, zs[ni], zp[ni], std::max(zd[ni], 0.3)};
        double ul2[4] = {0.0, zs[nj], zp[nj], std::max(zd[nj], 0.3)};
        int newk = std::min(a, b);
        int nk1 = newk + 1;
        for (int i = 1; i <= ia; ++i) {
            int iss = i;
            for (int j = 1; j <= b + 1; ++j) {
                int jss = j;
                for (int k = 1; k <= nk1; ++k) {
                    if (k > i || k > j) continue;
                    int kss = k;
                    int pi = std::max(npq[ni][i], iss);
                    int pj = std::max(npq[nj][j], jss);
                    S(s, i, j, k) =
                        ss_overlap(pi, pj, iss, jss, kss, ul1[i], ul2[j], r, a0);
                }
            }
        }
    }
    for (int i = 1; i <= ia; ++i) {  // i over L = s, p, d for atom a
        int kmin = 4 - i, kmax = 2 + i;
        for (int j = 1; j <= ib; ++j) {  // j over L = s, p, d for atom b
            double aa, bb;
            if (j == 2) {
                aa = -1.0;
                bb = 1.0;
            } else {
                aa = 1.0;
                if (j == 3)
                    bb = -1.0;
                else
                    bb = 1.0;
            }
            int lmin = 4 - j, lmax = 2 + j;
            for (int k = kmin; k <= kmax; ++k) {  // M_L slots for atom a
                for (int l = lmin; l <= lmax; ++l) {  // M_L slots for atom b
                    int ii = ival[i][k], jj = ival[j][l];
                    if (ii == 0 || jj == 0) continue;
                    di[ii][jj] = S(s, i, j, 1) * (cval(c.data(), i, k, 3) * cval(c.data(), j, l, 3)) * aa +
                                 S(s, i, j, 2) * (cval(c.data(), i, k, 4) * cval(c.data(), j, l, 4) +
                                                  cval(c.data(), i, k, 2) * cval(c.data(), j, l, 2)) * bb +
                                 S(s, i, j, 3) * (cval(c.data(), i, k, 5) * cval(c.data(), j, l, 5) +
                                                  cval(c.data(), i, k, 1) * cval(c.data(), j, l, 1));
                }
            }
        }
    }
    if (ni == 8 && nj == 1) {
        fprintf(stderr, "[DIAT] O-H xyz=%+.6f %+.6f %+.6f S21=%+.8e S31=%+.8e S41=%+.8e\n",
                x2, y2, z2, di[2][1], di[3][1], di[4][1]);
        fprintf(stderr, "[DIAT] S(s,2,1,1..3)=%+.8e %+.8e %+.8e cval214=%+.8e cval133=%+.8e\n",
                S(s, 2, 1, 1), S(s, 2, 1, 2), S(s, 2, 1, 3),
                cval(c.data(), 2, 4, 3), cval(c.data(), 1, 3, 3));
        fflush(stderr);
    }
}
