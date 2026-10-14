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

#ifdef MOPAC_USE_OPENCL
#include "opencl_accel.h"
// Optional GPU offload (opt-in, off by default).  The accelerator is
// initialised once; if no usable OpenCL device is found it stays inert and
// the CPU paths below run unchanged (bitwise reference behaviour).
static OpenCLAccel& ocl_accel() {
    static OpenCLAccel a;
    static bool tried = false;
    if (!tried) { tried = true; a.init(); }
    return a;
}
#endif

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
#ifdef MOPAC_USE_OPENCL
            if (ocl_accel().available()) {
                // Opt-in GPU path: same upper-triangle rank-2 update, different
                // reduction order (see opencl_accel.h for the bitwise caveat).
                // cT2/cT1 are row-major n x k copies built above.
                std::vector<double> xmat((size_t)norbs * (size_t)norbs, 0.0);
                if (ocl_accel().rank2_update(norbs, k2, cT2.data(), 2.0 * sign,
                                             k1, cT1.data(), frac * sign, xmat.data())) {
                    for (int j = 1; j <= norbs; ++j)
                        for (int i = 1; i <= j; ++i) {
                            double v = xmat[(size_t)(i - 1) * norbs + (j - 1)];
                            pp[pack_idx(i, j)] = (i == j) ? v + cst : v;
                        }
                    break;
                }
                // GPU failure: fall back to the CPU loop below.
            }
#endif
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
#ifdef MOPAC_USE_OPENCL
            if (ocl_accel().available()) {
                // Opt-in GPU path: upper-triangle DSYRK with a row-major copy of
                // c's first ndubl columns; reduction order differs (see header).
                std::vector<double> xmat((size_t)norbs * (size_t)norbs, 0.0);
                if (ndubl > 0) {
                    std::vector<double> cT((size_t)norbs * (size_t)ndubl);
                    for (int i = 0; i < norbs; ++i)
                        for (int k = 0; k < ndubl; ++k)
                            cT[(size_t)i * ndubl + k] = c[(size_t)k * norbs + i];
                    if (!ocl_accel().syrk_upper(norbs, ndubl, occ, cT.data(), xmat.data()))
                        goto syrk_cpu;  // GPU failure -> CPU fallback
                }
                for (int j = 1; j <= norbs; ++j)
                    for (int i = 1; i <= j; ++i)
                        pp[pack_idx(i, j)] = xmat[(size_t)(i - 1) * norbs + (j - 1)];
                break;
            }
syrk_cpu:
#endif
            std::vector<double> xmat(norbs * norbs, 0.0);
#if 1 // MKL dsyrk: match official density_for_GPU.F90 iopc=4/5 path
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
            } else {
                // Official F90: dsyrk(k=ndubl=0, beta=0) zeroes the full matrix,
                // then dtrttp copies it to pp.  The old `if (ndubl > 0)` guard
                // left pp untouched (e.g. UHF beta density for H/Li/Na single
                // atoms kept the initial 0.5 guess, corrupting the Fock build).
                int pN = (norbs * (norbs + 1)) / 2;
                for (int idx = 1; idx <= pN; ++idx) pp[idx] = 0.0;
            }
            int pN = (norbs * (norbs + 1)) / 2;
            break;
#endif
        }
        default:
            break;
    }
}
