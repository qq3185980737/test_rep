// rsp.cpp — symmetric eigensolver (packed lower triangle a(n*(n+1)/2)).
// Matches official rsp.F90: n==1 special case, otherwise forwards to
// eigenvectors_LAPACK.  NOTE (2026-10-06): forwarding to the default scalar
// Jacobi made the C12H26 optimisation diverge (2000 cycles, no convergence,
// -63.03), so this file keeps its own explicit MKL dsyevd('V','U') path —
// the EF Hessian diagonalisation is direction-sensitive and needs LAPACK's
// tridiagonal-reduction eigenvector phases; the scalar Jacobi (bit-exact for
// SCF single points) picks different near-degenerate eigenvector directions
// for EF Hessians (6 translation/rotation zero-modes), which corrupts the
// fx projection.  MOPAC_USE_MKL_DIAG=1 routes eigenvectors_LAPACK to dsyevd
// for the SCF path; rsp keeps its own dsyevd regardless.
#include "rsp.h"
#include <cmath>
#include <cstdio>
#include <vector>

#ifdef MOPAC_USE_MKL
#include <mkl_lapacke.h>
extern "C" void dsyevd_(char* jobz, char* uplo, int* n, double* a, int* lda,
                        double* w, double* work, int* lwork, int* iwork,
                        int* liwork, int* info);
#endif

void rsp(double* a, int n, double* root, double* vect) {
    if (n == 1) { root[0] = a[0]; vect[0] = 1.0; return; }
#ifdef MOPAC_USE_MKL
    // Unpack packed lower triangle into full symmetric matrix vect
    // (column-major).  Packed order (Fortran row-major lower):
    // (1,1),(2,1),(2,2),(3,1),(3,2),(3,3),...
    // Diagonal perturbation (i+1)*1e-10 (0-based i == Fortran i*1e-10).
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) vect[j * n + i] = 0.0;
    int k = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double v = a[k];
            if (i == j) v += (i + 1) * 1.0e-10;
            vect[j * n + i] = v;   // row i, col j (column-major storage)
            vect[i * n + j] = v;
            ++k;
        }
    std::vector<double> w(n);
    // Call Fortran dsyevd exactly like the reference eigenvectors_LAPACK.F90:
    // two calls (workspace query, then actual), jobz='v', uplo='u'.
    char jobz = 'v', uplo = 'u';
    int lwork = -1, liwork = -1, info = 0, lda = n;
    double work_q = 0.0;
    int iwork_q = 0;
    dsyevd_(&jobz, &uplo, &n, vect, &lda, w.data(), &work_q, &lwork,
            &iwork_q, &liwork, &info);
    if (info != 0) {
        for (int i = 0; i < n; ++i) root[i] = 0.0;
        for (int i = 0; i < n * n; ++i) vect[i] = 0.0;
        return;
    }
    lwork = (int)work_q;
    liwork = iwork_q;
    std::vector<double> work(lwork, 0.0);
    std::vector<int> iwork(liwork, 0);
    dsyevd_(&jobz, &uplo, &n, vect, &lda, w.data(), work.data(), &lwork,
            iwork.data(), &liwork, &info);
    if (info != 0) {
        for (int i = 0; i < n; ++i) root[i] = 0.0;
        for (int i = 0; i < n * n; ++i) vect[i] = 0.0;
        return;
    }
    for (int i = 0; i < n; ++i) root[i] = w[i];
#if 0  // CANON-DEG experiment (paper ①, CLOSED 2026-10-06): deterministic
    // degenerate-basis canonicalization.  Canonicalization itself is correct
    // (only rewrites degenerate columns; hessc/ec/grad unaffected bit-wise),
    // BUT fx[i] = u-row_i . grad contains components of the degenerate
    // eigenvectors on every DOF, so ANY canonicalization changes all fx and
    // the EF trajectory (C12H26 227 -> 1253 steps, worse convergence); with
    // DSPG grad-projection 918 steps.  Skipping degenerate columns or using
    // column-form fx diverges outright (G=72.9 / G=550).  Conclusion: EF's
    // row-semantics fx is intrinsically sensitive to the degenerate-basis
    // choice — a paper "discussion" finding, not a method.  Keep OFF for the
    // bit-level-official path.
    // CONCLUSION (2026-10-06): the canonicalization itself is correct (only
    // rewrites the degenerate columns; hessc/ec/non-degenerate columns and grad
    // all stay bit-identical), BUT fx[i] = u-row_i . grad contains components of
    // the degenerate eigenvectors on every DOF, so rewriting the degenerate
    // basis changes ALL fx[i>=7] and therefore the step d and the whole EF
    // trajectory (C12H26: 227 -> 1253 steps, worse convergence).  Kept off for
    // the bit-level-official path; re-enable only under the paper branch
    // (deterministic degenerate-basis method), ideally after projecting grad
    // onto the non-degenerate subspace so fx becomes insensitive to the
    // degenerate basis choice.
    // (see task log).  Disabled for the A/B rollback experiment; re-enable by
    // flipping to #if 1 and rebuilding.
    // --- Deterministic canonicalization of the near-degenerate subspace ---
    // Goal: make eigenvectors of degenerate (|lambda|<1e-6) modes depend only
    // on the subspace itself (i.e. on the input matrix), not on LAPACK's
    // internal basis choice.  Projector P = V_d V_d^T; pick canonical
    // directions by projecting standard basis vectors e_j in fixed order and
    // greedy orthonormalization.  Result is continuous in the input (no
    // basis-jump), so two inputs differing by 1 ULP give canonical bases
    // differing continuously, not by arbitrary rotation.
    {
        const double tol = 1.0e-6;
        std::vector<int> deg;
        for (int i = 0; i < n; ++i)
            if (std::fabs(w[i]) < tol) deg.push_back(i);
        int m = (int)deg.size();
        if (m > 0 && m < n && n > 50) {  // EF Hessian only; small calls
                                         // (molsymy 3x3, axis, etc.) must be
                                         // left untouched — canonicalizing their
                                         // degenerate vectors corrupts point
                                         // group / symmetry handling.
            std::fprintf(stderr, "[CANON] m=%d deg=", m);
            for (int k = 0; k < m; ++k) std::fprintf(stderr, " %d", deg[k]);
            std::fprintf(stderr, "\n w0-7:");
            for (int k = 0; k < 8 && k < n; ++k)
                std::fprintf(stderr, " %.6e", w[k]);
            std::fprintf(stderr, "\n");
            // Debug: col0/col1 before any rewrite (col-major: vect[r*n] etc.)
            std::fprintf(stderr, "[CANON] in0:");
            for (int r = 0; r < 4; ++r) std::fprintf(stderr, " %.9e", a[r]);
            std::fprintf(stderr, "\n[CANON] col0pre:");
            for (int r = 0; r < 4; ++r) std::fprintf(stderr, " %.9e", vect[r * n]);
            std::fprintf(stderr, "\n[CANON] col1pre:");
            for (int r = 0; r < 4; ++r)
                std::fprintf(stderr, " %.9e", vect[1 + r * n]);
            std::fprintf(stderr, "\n[CANON] col7pre:");
            for (int r = 0; r < 4; ++r)
                std::fprintf(stderr, " %.9e", vect[7 + r * n]);
            std::fprintf(stderr, "\n");
            std::vector<double> Q(n * m, 0.0);
            int qc = 0;
            for (int j = 0; j < n && qc < m; ++j) {
                std::vector<double> wj(n, 0.0);
                for (int k = 0; k < m; ++k) {
                    double c = vect[deg[k] + j * n];  // V_d^T e_j
                    for (int r = 0; r < n; ++r)
                        wj[r] += c * vect[deg[k] + r * n];
                }
                for (int q = 0; q < qc; ++q) {
                    double dot = 0.0;
                    for (int r = 0; r < n; ++r) dot += wj[r] * Q[q * n + r];
                    for (int r = 0; r < n; ++r) wj[r] -= dot * Q[q * n + r];
                }
                double nrm = 0.0;
                for (int r = 0; r < n; ++r) nrm += wj[r] * wj[r];
                nrm = std::sqrt(nrm);
                if (nrm > 1e-12) {
                    for (int r = 0; r < n; ++r) Q[qc * n + r] = wj[r] / nrm;
                    ++qc;
                }
            }
            if (qc == m) {
                for (int k = 0; k < m; ++k)
                    for (int r = 0; r < n; ++r)
                        vect[deg[k] + r * n] = Q[k * n + r];
                std::fprintf(stderr, "[CANON] col0post:");
                for (int r = 0; r < 4; ++r)
                    std::fprintf(stderr, " %.9e", vect[r * n]);
                std::fprintf(stderr, "\n[CANON] col1post:");
                for (int r = 0; r < 4; ++r)
                    std::fprintf(stderr, " %.9e", vect[1 + r * n]);
                std::fprintf(stderr, "\n[CANON] col7post:");
                for (int r = 0; r < 4; ++r)
                    std::fprintf(stderr, " %.9e", vect[7 + r * n]);
                std::fprintf(stderr, "\n");
            }
        }
    }
#endif  // CANON-DEG experiment
#else
    // No MKL: forward to eigenvectors_LAPACK (scalar Jacobi).  Only used when
    // building without Intel oneAPI; not bit-exact for EF Hessians.
    extern void eigenvectors_LAPACK(double* eigenvecs, double* xmat,
                                    double* eigvals, int ndim);
    std::vector<double> r(n + 1, 0.0);
    eigenvectors_LAPACK(vect, a, r.data(), n);
    for (int i = 0; i < n; ++i) root[i] = r[i + 1];
#endif  // MOPAC_USE_MKL
}

// 1-based Fortran-style wrapper for callers using std::vector (e.g. powsq.cpp).
void rsp(const std::vector<double>& a, int n, std::vector<double>& eig,
         std::vector<double>& pvec) {
    rsp(const_cast<double*>(a.data() + 1), n, eig.data() + 1, pvec.data() + 1);
}
