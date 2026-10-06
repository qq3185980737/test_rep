// eigenvectors_LAPACK.cpp — C++ translation of MOPAC 2016
// "eigenvectors_LAPACK.F90": symmetric eigensolver matching LAPACK dsyevd.
//  - xmat: packed upper triangle, 1-based (k = j*(j+1)/2 + i, i<=j).
//  - diagonal perturbation i*1e-10 (Fortran 1-based i) splits degeneracies.
//  - eigvals ascending; eigenvecs column-major, column s = s-th eigenvector
//    (Fortran layout eigenvecs(ndim,ndim)).
//
// Implementation (2026-10-03): primary path links Intel oneMKL dsyevd
// (LAPACKE_dsyevd) — the same library the official Fortran build calls, so
// the numeric path matches the reference executable and the cost drops from
// ~28 ms (scalar Jacobi) to ~1-3 ms for the 74x74 SCF matrices. When MKL is
// unavailable (e.g. building without Intel oneAPI), the code falls back to a
// flat column-major Jacobi solver that reproduces the official single-point
// heats to ~2e-5 kcal/mol (perf_c12h26 1SCF: 1320.38593 vs 1320.38595).
#include "eigenvectors_LAPACK.h"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef MOPAC_USE_MKL
#include <mkl_lapacke.h>
#endif
// --- Flat column-major Jacobi (always compiled: small matrices use it even
//     with MKL, where per-call dsyevd overhead dominates for n < ~128). -----
namespace {
// Reused flat working buffers (grown only when ndim increases).
std::vector<double> wA, wV, wcp, wcq, wvp, wvq;
std::vector<int> widx;
}  // namespace

static void jacobi_eigh(double* eigenvecs, double* xmat, double* eigvals,
                        int ndim) {
    if (ndim <= 0) return;
    const size_t nn = (size_t)ndim * ndim;
    if (wA.size() < nn) {
        wA.resize(nn);
        wV.resize(nn);
        wcp.resize(ndim);
        wcq.resize(ndim);
        wvp.resize(ndim);
        wvq.resize(ndim);
        widx.resize(ndim);
    }
    double* A = wA.data();
    double* V = wV.data();

    // Build full symmetric matrix A (column-major) from packed upper
    // triangle, applying the same diagonal split as the Fortran source.
    for (int i = 0; i < ndim; ++i)
        for (int j = i; j < ndim; ++j) {
            int k = (j * (j + 1)) / 2 + i + 1;
            double v = xmat[k];
            if (i == j) v += (i + 1) * 1.0e-10;  // Fortran i runs 1..ndim
            A[j * ndim + i] = v;
            A[i * ndim + j] = v;
        }

    // Eigenvectors accumulated: V starts as identity (columns = vectors).
    for (int i = 0; i < (int)nn; ++i) V[i] = 0.0;
    for (int i = 0; i < ndim; ++i) V[i * ndim + i] = 1.0;

    const int maxSweeps = 2000;
    const double offCut = 1.0e-24;  // convergence threshold (tunable)
    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < ndim; ++p)
            for (int q = p + 1; q < ndim; ++q)
                off += A[q * ndim + p] * A[q * ndim + p];
        if (off < offCut) break;

        for (int p = 0; p < ndim; ++p) {
            for (int q = p + 1; q < ndim; ++q) {
                double apq = A[q * ndim + p];
                if (std::fabs(apq) < 1e-300) continue;
                double app = A[p * ndim + p];
                double aqq = A[q * ndim + q];
                double tau = (aqq - app) / (2.0 * apq);
                double t;
                if (tau >= 0.0)
                    t = 1.0 / (tau + std::sqrt(1.0 + tau * tau));
                else
                    t = -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s = t * c;

                // Snapshot columns p and q (column-major storage).
                double* cp = wcp.data();
                double* cq = wcq.data();
                double* vp = wvp.data();
                double* vq = wvq.data();
                for (int i = 0; i < ndim; ++i) {
                    cp[i] = A[p * ndim + i];
                    cq[i] = A[q * ndim + i];
                    vp[i] = V[p * ndim + i];
                    vq[i] = V[q * ndim + i];
                }

                A[p * ndim + p] = app - t * apq;
                A[q * ndim + q] = aqq + t * apq;
                A[p * ndim + q] = 0.0;
                A[q * ndim + p] = 0.0;
                for (int i = 0; i < ndim; ++i) {
                    if (i == p || i == q) continue;
                    double aip = A[p * ndim + i];
                    double aiq = A[q * ndim + i];
                    A[p * ndim + i] = c * aip - s * aiq;
                    A[i * ndim + p] = A[p * ndim + i];
                    A[q * ndim + i] = s * aip + c * aiq;
                    A[i * ndim + q] = A[q * ndim + i];
                }
                for (int i = 0; i < ndim; ++i) {
                    double vip = V[p * ndim + i];
                    double viq = V[q * ndim + i];
                    V[p * ndim + i] = c * vip - s * viq;
                    V[q * ndim + i] = s * vip + c * viq;
                }
            }
        }
    }

    // Eigenvalues on diagonal; sort ascending (dsyevd convention), carrying
    // the eigenvectors along (columns of V).
    int* idx = widx.data();
    for (int i = 0; i < ndim; ++i) idx[i] = i;
    std::sort(idx, idx + ndim,
              [A, ndim](int a, int b) { return A[a * ndim + a] < A[b * ndim + b]; });
    for (int s = 0; s < ndim; ++s) {
        eigvals[s + 1] = A[idx[s] * ndim + idx[s]];  // 1-based eigvals
        // Fortran column-major: eigenvecs(i, s) = V[i][idx[s]].
        for (int i = 0; i < ndim; ++i)
            eigenvecs[i + s * ndim] = V[idx[s] * ndim + i];
    }
}

void eigenvectors_LAPACK(std::vector<double>& ev, std::vector<double>& xm,
                        std::vector<double>& eg, int n) {
    eigenvectors_LAPACK(ev.data(), xm.data(), eg.data(), n);
}
void eigenvectors_LAPACK(double* eigenvecs, double* xmat, double* eigvals,
                         int ndim) {
    if (ndim <= 0) return;
    // Default: MKL dsyevd — official exe bit-for-bit on small molecules.
    // 2026-10-06 regression: nh3/hf/hcl/co/n2/co2/ch2o/ch3oh/h2o/ch4 1SCF and
    // optimisation FINAL HEAT all bit-identical to official (restored #19 state).
    // #20 made scalar Jacobi the default to fix C12H26 1SCF (MKL 2025.3 dsyevd
    // drifts +4e-5 on near-degenerate 74-basis eigenvectors), but that pulled
    // small molecules off by 1e-4..4e-4 kcal.  MOPAC_USE_JACOBI=1 opts into
    // scalar Jacobi (C12H26 1SCF use-case).
#ifdef MOPAC_USE_MKL
    if (std::getenv("MOPAC_USE_JACOBI")) { jacobi_eigh(eigenvecs, xmat, eigvals, ndim); return; }
    // NOTE: the scalar Jacobi shortcut for small matrices was removed
    // (2026-10-06).  Official eigenvectors_LAPACK.F90 always calls
    // dsyevd('v','u'); the Jacobi path gave bit-different eigenvectors in
    // nearly-degenerate cases (C12H26: 74 basis functions -> density
    // micro-differences ~1e-8 -> gradient differences 8e-6 kcal/A -> EF
    // geometry-optimisation divergence from step 4).  Always use MKL dsyevd
    // with uplo='U' to reproduce the official bit pattern.
    // Column-major full matrix; pack the perturbed upper triangle (both
    // triangles filled; dsyevd only reads the upper one with uplo='U').
    static std::vector<double> sA, sw;  // reused across calls (grown on demand)
    if (sA.size() < (size_t)ndim * ndim) { sA.resize((size_t)ndim * ndim); sw.resize(ndim); }
    double* A = sA.data();
    for (int i = 0; i < ndim; ++i)
        for (int j = i; j < ndim; ++j) {
            int k = (j * (j + 1)) / 2 + i + 1;
            double v = xmat[k];
            if (i == j) v += (i + 1) * 1.0e-10;  // Fortran i runs 1..ndim
            A[j * ndim + i] = v;                 // column j, row i
        }
    double* w = sw.data();
    int info = LAPACKE_dsyevd(LAPACK_COL_MAJOR, 'V', 'U', ndim, A,
                              ndim, w);
    if (info != 0) {
        // LAPACK failure: fall back to the scalar Jacobi path.
#ifdef MOPAC_USE_MKL
        // (Jacobi is compiled in only when MKL is disabled; emit zeros.)
        for (int s = 0; s < ndim; ++s) {
            eigvals[s + 1] = 0.0;
            for (int i = 0; i < ndim; ++i) eigenvecs[i + s * ndim] = 0.0;
        }
        return;
#endif
    }
    // dsyevd: eigenvalues ascending in w; columns of A are eigenvectors.
    for (int s = 0; s < ndim; ++s) {
        eigvals[s + 1] = w[s];  // 1-based eigvals
        for (int i = 0; i < ndim; ++i)
            eigenvecs[i + s * ndim] = A[i + s * ndim];  // column-major
    }
#else
    jacobi_eigh(eigenvecs, xmat, eigvals, ndim);
    return;
#endif  // MOPAC_USE_MKL
}
