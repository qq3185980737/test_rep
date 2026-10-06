// diag_for_GPU.cpp — C++ translation of MOPAC 2016 "diag_for_GPU.F90".
// Pseudo-diagonalisation over the occupied/virtual MO intersection (CPU path;
// GPU branches execute the equivalent CPU computation).
#include "diag_for_GPU.h"
#include "molkst_C.h"
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef MOPAC_USE_MKL
#include <mkl.h>
#endif
#include <chrono>
#include <cstdio>

namespace molkst_C {
extern int num_threads;
}

namespace {
// idamax: index (0-based) of element with max |x|.
int idamax(int n, const double* x, int incx) {
    int imax = 0;
    double vmax = -1.0;
    for (int i = 0, k = 0; i < n; ++i, k += incx) {
        double a = std::fabs(x[k]);
        if (a > vmax) { vmax = a; imax = i; }
    }
    return imax;
}
}  // namespace

void diag_for_GPU(std::vector<double>& fao, std::vector<double>& vector, int nocc,
                  std::vector<double>& eig, int mdim, int icalcn) {
    diag_for_GPU(fao.data(), vector.data(), nocc, eig.data(), mdim, icalcn);
}
void diag_for_GPU(double* fao, double* vector, int nocc, double* eig,
                  int norbs, int mpack) {
    if (nocc == norbs || nocc == 0) return;
#ifdef MOPAC_USE_MKL
    // Mirror official diag_for_GPU.F90 case(1): BLAS single-threaded.
    // Official sets mkl_set_num_threads(1) before the rotation loop; force
    // single-thread dgemm too for a clean bit-level experiment.
    mkl_set_num_threads(1);
#endif
    int n = norbs, mdim = norbs, lumo = nocc + 1, nvirt = n - nocc;
    int norbs2 = norbs * norbs;
    static long long _tu = 0, _tg = 0, _tt = 0, _tr = 0;
    static int _tc = 0;
    auto _a = std::chrono::steady_clock::now();
    static std::vector<double> fmo, fck;  // reused across calls (full overwrite below)
    const size_t need = (size_t)norbs * norbs;
    if (fmo.size() < need) { fmo.resize(need); fck.resize(need); }
    // Unpack fao (packed triangular; location by upper/lower rule).
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= n; ++j) {
            int in = (j > i) ? ((j * (j - 1)) / 2) + i : ((i * (i - 1)) / 2) + j;
            int ij = j + n * (i - 1);  // fmo(i,j) column-major, 1-based
            fmo[ij - 1] = fao[in];      // fao is 1-based (fao(1..mpack)); C++ f[0] unused
        }
    _tu += std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - _a).count();
    auto _b = std::chrono::steady_clock::now();
    // C(m,nn) = A(m,kk) * B(kk,nn); all column-major.
    // Cc(m,nn) = A^T(m,kk) * B(kk,nn); A supplied as A^T stored (m,kk) column-major,
    // i.e. A^T(ii,t) = A[ii*kk + t]. (F90 passes A=vector(lumo:n,:) as (kk,m)
    // column-major and uses 'T'; here we pass the transposed view instead.)
#ifdef MOPAC_USE_MKL
    // Serial reads A[ii*kk+t] for gemm_tn, which equals the column-major address
    // t + ii*kk of a (kk x m) matrix; so COLUMN-major 'T' matches the serial code
    // exactly (both are the same linear address).
    // NOTE: official diag_for_GPU.F90 calls dgemm unconditionally.  The serial
    // fallback below (< 500000 flops) gave a different summation order than MKL
    // dgemm, so fmo differed by a few ulps for ~200k-flop products (C12H26:
    // 74x74x37) and the pseudo-diagonalisation rotation set flipped at threshold
    // edges -> density micro-differences -> gradient differences (8e-6 kcal/A)
    // -> EF geometry-optimisation divergence from step 4 onward.  Always call
    // MKL dgemm to match the official bit pattern.
    auto gemm_nn = [&](int m, int kk, int nn, const double* A, const double* B, double* Cc) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, nn, kk,
                    1.0, A, m, B, kk, 0.0, Cc, m);
    };
    auto gemm_tn = [&](int m, int kk, int nn, const double* A, const double* B, double* Cc) {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, nn, kk,
                    1.0, A, kk, B, kk, 0.0, Cc, m);
    };
#else
    auto gemm_nn = [&](int m, int kk, int nn, const double* A, const double* B, double* Cc) {
        for (int jj = 0; jj < nn; ++jj)
            for (int ii = 0; ii < m; ++ii) {
                double acc = 0.0;
                for (int t = 0; t < kk; ++t) acc += A[t * m + ii] * B[jj * kk + t];
                Cc[jj * m + ii] = acc;
            }
    };
    auto gemm_tn = [&](int m, int kk, int nn, const double* A, const double* B, double* Cc) {
        for (int jj = 0; jj < nn; ++jj)
            for (int ii = 0; ii < m; ++ii) {
                double acc = 0.0;
                for (int t = 0; t < kk; ++t) acc += A[ii * kk + t] * B[jj * kk + t];
                Cc[jj * m + ii] = acc;
            }
    };
#endif
    if (nocc < nvirt) {
        // fck(n,nocc) = fmo(n,n) * vector(n,nocc)
        gemm_nn(n, n, nocc, fmo.data(), vector, fck.data());
        // fmo(nvirt,nocc) = vector(lumo:n,:)^T * fck
        static std::vector<double> vt;
        if (vt.size() < (size_t)n * nvirt) vt.resize((size_t)n * nvirt);
        auto _c = std::chrono::steady_clock::now();
        for (int c = 0; c < nvirt; ++c)
            for (int r = 0; r < n; ++r) vt[c * n + r] = vector[(lumo - 1 + c) * mdim + r];
        _tt += std::chrono::duration_cast<std::chrono::microseconds>(
                   std::chrono::steady_clock::now() - _c).count();
        gemm_tn(nvirt, n, nocc, vt.data(), fck.data(), fmo.data());
    } else {
        // fck(n,nvirt) = fmo(n,n) * vector(n, lumo:n)
        static std::vector<double> vs;
        if (vs.size() < (size_t)n * nvirt) vs.resize((size_t)n * nvirt);
        auto _c = std::chrono::steady_clock::now();
        for (int c = 0; c < nvirt; ++c)
            for (int r = 0; r < n; ++r) vs[c * n + r] = vector[(lumo - 1 + c) * mdim + r];
        _tt += std::chrono::duration_cast<std::chrono::microseconds>(
                   std::chrono::steady_clock::now() - _c).count();
        gemm_nn(n, n, nvirt, fmo.data(), vs.data(), fck.data());
        // fmo(nvirt,nocc) = fck^T * vector(n,nocc); fck is (n,nvirt) column-major
        // with column stride n (gemm_nn output), i.e. fck(row=t, col=ii) = fck[t + ii*n].
        // Build its transpose view (nvirt,n) column-major: fckT(ii,t) = fck(t,ii).
        static std::vector<double> fckT;
        if (fckT.size() < (size_t)nvirt * n) fckT.resize((size_t)nvirt * n);
        auto _c2 = std::chrono::steady_clock::now();
        for (int c = 0; c < nvirt; ++c)
            for (int r = 0; r < n; ++r) fckT[c * n + r] = fck[r + c * n];
        _tt += std::chrono::duration_cast<std::chrono::microseconds>(
                   std::chrono::steady_clock::now() - _c2).count();
        gemm_tn(nvirt, n, nocc, fckT.data(), vector, fmo.data());
    }
    _tg += std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - _b).count();
    int i = idamax(nocc * nvirt, fmo.data(), 1);
    double tiny = std::fabs(fmo[i]) * 0.05;
    double bigeps = 1.5e-7;
    int ij = 0;
    auto _d = std::chrono::steady_clock::now();
    for (int ii = 1; ii <= nocc; ++ii)
        for (int jj = lumo; jj <= n; ++jj) {
            double x = fmo[ij];
            ++ij;
            if (std::fabs(x) < tiny) continue;
            double a = eig[ii], b = eig[jj];
            double d = a - b;
            if (std::fabs(x / d) < bigeps) continue;
            double e = (d >= 0 ? 1.0 : -1.0) * std::sqrt(4.0 * x * x + d * d);
            double alpha = std::sqrt(0.5 * (1.0 + d / e));
            double beta = -(x >= 0 ? 1.0 : -1.0) * std::sqrt(1.0 - alpha * alpha);
            // Official: call drot(n, vector(1:n,ii), 1, vector(1:n,jj), 1, alpha, beta)
#ifdef MOPAC_USE_MKL
            cblas_drot(n, vector + (ii - 1) * mdim, 1, vector + (jj - 1) * mdim, 1, alpha, beta);
#else
            for (int k = 0; k < n; ++k) {
                double xi = vector[(ii - 1) * mdim + k];
                double yj = vector[(jj - 1) * mdim + k];
                vector[(ii - 1) * mdim + k] = alpha * xi + beta * yj;
                vector[(jj - 1) * mdim + k] = -beta * xi + alpha * yj;
            }
#endif
        }
    _tr += std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - _d).count();
    if (++_tc % 8 == 0)
        fprintf(stderr, "[DGFPROF] unpack=%.1f gemm=%.1f trans=%.1f rot=%.1fms n=%d\n",
                _tu / 1e3, _tg / 1e3, _tt / 1e3, _tr / 1e3, norbs);
#ifdef MOPAC_USE_MKL
    mkl_set_num_threads(molkst_C::num_threads);  // restore (official F90 L207)
#endif
}
