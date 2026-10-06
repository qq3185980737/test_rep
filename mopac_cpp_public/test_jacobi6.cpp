// test_jacobi6.cpp — verify Jacobi eigensolver residual per column.
#include "eigenvectors_LAPACK.h"
#include <cstdio>
#include <cmath>
#include <vector>

int main() {
    const int n = 6;
    double Mt[6][6] = {
        {1.0, 0.2, 0.0, 0.0, 0.0, 0.0},
        {0.2, 1.0, 0.1, 0.0, 0.0, 0.0},
        {0.0, 0.1, 1.0, 0.3, 0.0, 0.0},
        {0.0, 0.0, 0.3, 1.0, 0.1, 0.0},
        {0.0, 0.0, 0.0, 0.1, 1.0, 0.2},
        {0.0, 0.0, 0.0, 0.0, 0.2, 1.0}
    };
    std::vector<double> A(n * n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) {
            double s = 0;
            for (int k = 0; k < n; ++k) s += Mt[i][k] * Mt[j][k];
            A[i * n + j] = s;
        }
    std::vector<double> xmat(22, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            xmat[(j * (j + 1)) / 2 + i + 1] = A[i * n + j];
    std::vector<double> V(n * n, 0.0), lam(n, 0.0);
    eigenvectors_LAPACK(V.data(), xmat.data(), lam.data(), n);
    for (int c = 0; c < n; ++c) {
        double res = 0;
        for (int r = 0; r < n; ++r) {
            double av = 0;
            for (int k = 0; k < n; ++k) av += A[r * n + k] * V[k * n + c];
            double vd = V[r * n + c] * lam[c];
            res += (av - vd) * (av - vd);
        }
        std::printf("col %d lam=%.6f res=%.3e vec=[", c, lam[c]);
        for (int r = 0; r < n; ++r) std::printf(" %.4f", V[r * n + c]);
        std::printf(" ]\n");
    }
    return 0;
}
