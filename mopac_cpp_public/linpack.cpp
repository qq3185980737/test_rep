// linpack.cpp — LINPACK dgefa / dgedi, column-major double arrays (lda x n).
#include "linpack.h"

#include <cmath>

// Internal access: row i, col j (both 1-based) in column-major storage.
static inline double& A(double* a, int lda, int i, int j) {
    return a[(i - 1) + (j - 1) * lda];
}

void dgefa(double* a, int lda, int n, int* ipvt, int& info) {
    info = 0;
    int nm1 = n - 1;
    for (int k = 1; k <= nm1; ++k) {
        int kp1 = k + 1;
        // idamax over a(k:n, k): pick largest abs in column k rows k..n.
        int l = k;
        double best = -1.0;
        for (int i = k; i <= n; ++i) {
            double v = std::fabs(A(a, lda, i, k));
            if (v > best) { best = v; l = i; }
        }
        ipvt[k] = l;
        if (A(a, lda, l, k) == 0.0) {
            info = k;
        } else {
            if (l != k) {
                double t = A(a, lda, l, k);
                A(a, lda, l, k) = A(a, lda, k, k);
                A(a, lda, k, k) = t;
            }
            double t = -1.0 / A(a, lda, k, k);
            // dscal(n-k, t, &a(k+1,k))
            for (int i = k + 1; i <= n; ++i) A(a, lda, i, k) *= t;
            // row elimination
            for (int j = kp1; j <= n; ++j) {
                double t2 = A(a, lda, l, j);
                if (l != k) {
                    A(a, lda, l, j) = A(a, lda, k, j);
                    A(a, lda, k, j) = t2;
                }
                // daxpy(n-k, t2, a(k+1,k), 1, a(k+1,j), 1)
                // MKL daxpy uses FMA (single rounding); emulate with std::fma.
                for (int i = k + 1; i <= n; ++i)
                    A(a, lda, i, j) = std::fma(t2, A(a, lda, i, k),
                                               A(a, lda, i, j));
            }
        }
    }
    ipvt[n] = n;
    if (A(a, lda, n, n) == 0.0) info = n;
}

void dgedi(double* a, int lda, int n, const int* ipvt, double* det,
           double* work, int job) {
    // determinant
    if (job / 10 != 0) {
        det[0] = 1.0;   // det(1)
        det[1] = 0.0;   // det(2)
        const double ten = 10.0;
        for (int i = 1; i <= n; ++i) {
            if (ipvt[i] != i) det[0] = -det[0];
            det[0] = A(a, lda, i, i) * det[0];
            if (det[0] == 0.0) break;
            while (std::fabs(det[0]) < 1.0) { det[0] *= ten; det[1] -= 1.0; }
            while (std::fabs(det[0]) >= ten) { det[0] /= ten; det[1] += 1.0; }
        }
    }
    // inverse(U)
    // F90: if (Mod(job,10)==0) return  =>  inverse runs when 个位 != 0
    if (job % 10 == 0) return;
    for (int k = 1; k <= n; ++k) {
        A(a, lda, k, k) = 1.0 / A(a, lda, k, k);
        double t = -A(a, lda, k, k);
        // dscal(k-1, t, &a(1,k))
        for (int i = 1; i <= k - 1; ++i) A(a, lda, i, k) *= t;
        int kp1 = k + 1;
        for (int j = kp1; j <= n; ++j) {
            double tj = A(a, lda, k, j);
            A(a, lda, k, j) = 0.0;
            // daxpy(k, tj, a(1,k), 1, a(1,j), 1)
            for (int i = 1; i <= k; ++i)
                A(a, lda, i, j) = std::fma(tj, A(a, lda, i, k),
                                           A(a, lda, i, j));
        }
    }
    // inverse(U) * inverse(L)
    int nm1 = n - 1;
    if (nm1 < 1) return;
    for (int kb = 1; kb <= nm1; ++kb) {
        int k = n - kb;
        int kp1 = k + 1;
        for (int i = kp1; i <= n; ++i) {
            work[i] = A(a, lda, i, k);
            A(a, lda, i, k) = 0.0;
        }
        for (int j = kp1; j <= n; ++j) {
            double tj = work[j];
            // daxpy(n, tj, a(1,j), 1, a(1,k), 1)
            for (int i = 1; i <= n; ++i)
                A(a, lda, i, k) = std::fma(tj, A(a, lda, i, j),
                                           A(a, lda, i, k));
        }
        int l = ipvt[k];
        if (l != k) {
            // dswap(n, a(1,k), 1, a(1,l), 1)
            for (int i = 1; i <= n; ++i) {
                double t = A(a, lda, i, k);
                A(a, lda, i, k) = A(a, lda, i, l);
                A(a, lda, i, l) = t;
            }
        }
    }
}
