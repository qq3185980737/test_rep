// mult_symm_AB.cpp — C++ translation of MOPAC 2016 "mult_symm_AB.F90".
// Symmetric-packed matrix multiply C = alpha*A*B + beta*C (packed storage).
//   iopc=1: mamult (packed in-place); iopc=3: unpack->dgemm->repack (CPU).
//   iopc=2/4: GPU paths — executed here with the equivalent CPU dgemm path.
#include "mult_symm_AB.h"
#include "mamult.h"
#include <cmath>
#include <vector>

extern void mamult(const double* a, const double* b, double* c, int n, double one);

namespace {
// Unpack upper triangle of a packed (m*(m+1)/2) symmetric matrix into full m x m,
// column-major; mirror the upper part into the lower triangle.
void dtpttr_full(int n, const double* packed, std::vector<double>& full) {
    int k = 0;
    for (int j = 1; j <= n; ++j) {
        for (int i = 1; i <= j; ++i) {
            full[(j - 1) * n + (i - 1)] = packed[k];
            ++k;
        }
    }
    for (int i = 1; i <= n - 1; ++i)
        for (int j = i + 1; j <= n; ++j)
            full[(i - 1) * n + (j - 1)] = full[(j - 1) * n + (i - 1)];
}
// Repack upper triangle of full symmetric matrix into packed (column-major upper).
void dtrttp_pack(int n, const std::vector<double>& full, double* packed) {
    int k = 0;
    for (int j = 1; j <= n; ++j)
        for (int i = 1; i <= j; ++i) {
            packed[k] = full[(j - 1) * n + (i - 1)];
            ++k;
        }
}
}  // namespace

void mult_symm_AB(double* a, double* b, double alpha, int ndim, int mdim,
                  double* c, double beta, int iopc) {
    if (iopc == 1) {
        mamult(a, b, c, ndim, beta);
        return;
    }
    int n = ndim;
    std::vector<double> xa(n * n, 0.0), xb(n * n, 0.0), xc(n * n, 0.0);
    dtpttr_full(n, a, xa);
    dtpttr_full(n, b, xb);
    if (beta != 0.0) dtpttr_full(n, c, xc);
    // dgemm("N","N",n,n,n,alpha,xa,n,xb,n,beta,xc,n) — column-major A*B.
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            double acc = beta * xc[j * n + i];
            for (int k = 0; k < n; ++k) acc += alpha * xa[k * n + i] * xb[j * n + k];
            xc[j * n + i] = acc;
        }
    dtrttp_pack(n, xc, c);
}
