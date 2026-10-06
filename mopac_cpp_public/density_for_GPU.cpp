// density_for_GPU.cpp — C++ translation of MOPAC 2016 "density_for_GPU.F90".
// Builds the (packed) density matrix P from the occupied eigenvectors:
//   P = 2*c_occ*c_occ^T (+ fractional part), or the positron equivalent,
// then packs the upper triangle (dtrttp) into pp.
// The GPU-only iopc=2/4 paths (CUBLAS) are implemented with the same CPU
// arithmetic so behaviour matches the BLAS paths (iopc=3/5).
#include "density_for_GPU.h"

#include <cmath>
#include <cstdio>
#include <vector>

#ifdef MOPAC_USE_MKL
#include <mkl.h>
#endif
#include <omp.h>

namespace {
inline int pack_idx(int i, int j) {  // upper triangle, Fortran packed (1-based i<=j)
    return j * (j - 1) / 2 + i;      // pp(pack_idx(i,j)) is P(i,j); caller uses pp[pack_idx-1]
}
}

void density_for_GPU(std::vector<double>& c, double fract, int ndubl, int nsingl,
                     double en, int nmos, int iopen, int nclose,
                     std::vector<double>& p, int icalcn) {
    density_for_GPU(c.data(), fract, ndubl, nsingl, en, nmos, iopen, nclose, p.data(), icalcn);
}
void density_for_GPU(double* c, double fract, int ndubl, int nsingl,
                     double occ, int /*mpack*/, int norbs, int mode, double* pp,
                     int iopc) {
    double sign, frac, cst;
    int nl1, nl2, nu1, nu2;
    if (ndubl != 0 && nsingl > (norbs / 2) && mode != 2) {
        // Positron equivalent.
        sign = -1.0;
        frac = occ - fract;
        cst = occ;
        nl2 = nsingl + 1;
        nu2 = norbs;
        nl1 = ndubl + 1;
        nu1 = nsingl;
    } else {
        // Electron equivalent.
        sign = 1.0;
        frac = fract;
        cst = 0.0;
        nl2 = 1;
        nu2 = ndubl;
        nl1 = ndubl + 1;
        nu1 = nsingl;
    }
    switch (iopc) {
        case 2:
        case 3: {  // gemm path: xmat = 2*sign*c2*c2^T + frac*sign*c1*c1^T; diag += cst
            int nl21 = nl2 < norbs ? nl2 : norbs;
            int nl11 = nl1 < norbs ? nl1 : norbs;
            // Serial semantics:  xmat(i,j) = 2*sign*Sum_k2 c(k,i)c(k,j)
            //                              + frac*sign*Sum_k1 c(k,i)c(k,j)  (+cst on diag)
            // Only the upper triangle is needed (packed below), so compute just
            // i<=j; the two rank-k sums are accumulated in the same arithmetic
            // order as the serial code, keeping results bit-identical.  The
            // transposed cache cT[i][k]=c[nl+k][i] makes the inner k loop read
            // contiguous memory.  Each (i,j) is completed by exactly one thread.
            int k2 = nu2 - nl21 + 1, k1 = nu1 - nl11 + 1;
            static std::vector<double> cT2, cT1;
            const size_t cap2 = (size_t)norbs * (size_t)std::max(k2, 1);
            const size_t cap1 = (size_t)norbs * (size_t)std::max(k1, 1);
            if (cT2.size() < cap2) cT2.resize(cap2);
            if (cT1.size() < cap1) cT1.resize(cap1);
            if (k2 > 0) {
#pragma omp parallel for schedule(static)
                for (int i = 1; i <= norbs; ++i) {
                    double* ci = &cT2[(size_t)(i - 1) * k2];
                    for (int k = 0; k < k2; ++k)
                        ci[k] = c[(size_t)(nl21 - 1 + k) * norbs + (i - 1)];
                }
            }
            if (k1 > 0) {
#pragma omp parallel for schedule(static)
                for (int i = 1; i <= norbs; ++i) {
                    double* ci = &cT1[(size_t)(i - 1) * k1];
                    for (int k = 0; k < k1; ++k)
                        ci[k] = c[(size_t)(nl11 - 1 + k) * norbs + (i - 1)];
                }
            }
#pragma omp parallel for schedule(static) collapse(2)
            for (int j = 1; j <= norbs; ++j)
                for (int i = 1; i <= j; ++i) {
                    double s = 0.0;
                    if (k2 > 0) {
                        const double* ci = &cT2[(size_t)(i - 1) * k2];
                        const double* cj_ = &cT2[(size_t)(j - 1) * k2];
                        for (int k = 0; k < k2; ++k) s += ci[k] * cj_[k];
                        s *= 2.0 * sign;
                    }
                    if (k1 > 0) {
                        const double* ci = &cT1[(size_t)(i - 1) * k1];
                        const double* cj_ = &cT1[(size_t)(j - 1) * k1];
                        double s1 = 0.0;
                        for (int k = 0; k < k1; ++k) s1 += ci[k] * cj_[k];
                        s += s1 * (frac * sign);  // serial order: xmat += (frac*sign)*s1
                    }
                    pp[pack_idx(i, j)] = (i == j) ? s + cst : s;
                }
            break;
        }
        case 4:
        case 5: {  // syrk path: xmat = occ*c(:,1:ndubl)*c(:,1:ndubl)^T
            std::vector<double> xmat(norbs * norbs, 0.0);
#if 0 // MKL temporarily off for bisection
            // c is column-major (LAPACK output; column = MO).  Serial code fills
            // xmat in C-style row-major linearisation ((i,j)->(i-1)*n+(j-1)).
            // We call dsyrk COLUMN-major (reads the same K columns of c as the
            // serial loop) but its result lands at the transposed location; P is
            // symmetric, so packing reads the transposed index instead.
            if (ndubl > 0)
                cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, norbs, ndubl,
                            occ, c, norbs, 1.0, xmat.data(), norbs);
            for (int j = 1; j <= norbs; ++j)
                for (int i = 1; i <= j; ++i)
                    pp[pack_idx(i, j)] = xmat[(size_t)(j - 1) * norbs + (i - 1)];
            {   // DEBUG: C12/C80 dsyrk check
                static int dbg = 0;
                if (dbg++ < 2) {
                    std::vector<double> xm2(norbs * norbs, 0.0);
                    for (int i = 1; i <= norbs; ++i)
                        for (int j = 1; j <= norbs; ++j) {
                            double s = 0.0;
                            for (int k = 1; k <= ndubl; ++k)
                                s += c[(k - 1) * norbs + (i - 1)] * c[(k - 1) * norbs + (j - 1)];
                            xm2[(i - 1) * norbs + (j - 1)] = occ * s;
                        }
                    double mx = 0.0; int mi = -1, mj = -1;
                    for (int j = 1; j <= norbs; ++j)
                        for (int i = 1; i <= j; ++i) {
                            double d = std::fabs(xm2[(i - 1) * norbs + (j - 1)] - xmat[(j - 1) * norbs + (i - 1)]);
                            if (d > mx) { mx = d; mi = i; mj = j; }
                        }
                    // compare PACKED pp vs serial-packed
                    double mp = 0.0; int mip = -1;
                    for (int j = 1; j <= norbs; ++j)
                        for (int i = 1; i <= j; ++i) {
                            int idx = pack_idx(i, j);
                            double d = std::fabs(xm2[(i - 1) * norbs + (j - 1)] - pp[idx]);
                            if (d > mp) { mp = d; mip = idx; }
                        }
                    fprintf(stderr, "[DBG] dsyrk-check maxdiff=%.3e at(%d,%d) pp_maxdiff=%.3e at(%d) norbs=%d k=%d occ=%.1f\n",
                            mx, mi, mj, mp, mip, norbs, ndubl, occ);
                    fflush(stderr);
                }
            }
            int pN = (norbs * (norbs + 1)) / 2;
            break;
#else
            // Transposed cache cT[i][k] = c[k][i] makes the inner k loop read
            // contiguous memory.  Each (i,j) pair is accumulated by exactly one
            // thread in the identical k order as the serial code, so the result
            // is bit-identical; only the upper triangle (i<=j) is computed.
            static std::vector<double> cT;
            const size_t cTcap = (size_t)norbs * (size_t)std::max(ndubl, 1);
            if (cT.size() < cTcap) cT.resize(cTcap);
            if (ndubl > 0) {
#pragma omp parallel for schedule(static)
                for (int i = 1; i <= norbs; ++i) {
                    double* ci = &cT[(size_t)(i - 1) * ndubl];
                    for (int k = 0; k < ndubl; ++k)
                        ci[k] = c[(size_t)k * norbs + (i - 1)];
                }
#pragma omp parallel for schedule(static) collapse(2)
                for (int j = 1; j <= norbs; ++j)
                    for (int i = 1; i <= j; ++i) {
                        double s = 0.0;
                        const double* ci = &cT[(size_t)(i - 1) * ndubl];
                        const double* cj_ = &cT[(size_t)(j - 1) * ndubl];
                        for (int k = 0; k < ndubl; ++k) s += ci[k] * cj_[k];
                        pp[pack_idx(i, j)] = occ * s;
                    }
            }
            int pN = (norbs * (norbs + 1)) / 2;
            break;
#endif
        }
        default:
            break;
    }
}
