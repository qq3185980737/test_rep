// test_density_for_GPU.cpp — verify density_for_GPU against hand-computed
// values for norbs=2 with identity eigenvectors.
#include "density_for_GPU.h"
#include <cstdio>
#include <cmath>

static bool near(double a, double b, double tol = 1e-12) {
    return std::fabs(a - b) <= tol;
}

int main() {
    bool ok = true;
    // c = I (column-major): c[0]=1 (row1,col1), c[3]=1 (row2,col2)
    double c[4] = {1.0, 0.0, 0.0, 1.0};
    double pp[3];
    int norbs = 2, nlower = 3;

    // --- Electron equivalent, iopc=3 (dgemm path), closed shell ---
    // ndubl=1, nsingl=0: P = 2*c1*c1^T = diag(2, 0) -> packed (2, 0, 0)
    for (int i = 0; i < 3; ++i) pp[i] = -9.0;
    density_for_GPU(c, 1.0, 1, 0, 2.0, nlower, norbs, 0, pp, 3);
    if (!(near(pp[0], 2.0) && near(pp[1], 0.0) && near(pp[2], 0.0))) {
        std::printf("FAIL electron dgemm: %g %g %g\n", pp[0], pp[1], pp[2]);
        ok = false;
    }

    // --- Positron equivalent, iopc=3 ---
    // ndubl=1, nsingl=2 (>norbs/2): sign=-1, frac=occ-fract=1, cst=2
    // xmat = -2*c2*c2^T - 1*c2*c2^T + 2*I = diag(2, -1) -> packed (2, 0, -1)
    for (int i = 0; i < 3; ++i) pp[i] = -9.0;
    density_for_GPU(c, 1.0, 1, 2, 2.0, nlower, norbs, 0, pp, 3);
    if (!(near(pp[0], 2.0) && near(pp[1], 0.0) && near(pp[2], -1.0))) {
        std::printf("FAIL positron dgemm: %g %g %g\n", pp[0], pp[1], pp[2]);
        ok = false;
    }

    // --- iopc=5 (dsyrk path): P = occ*c(:,1:ndubl)*c(:,1:ndubl)^T = diag(2,0) ---
    for (int i = 0; i < 3; ++i) pp[i] = -9.0;
    density_for_GPU(c, 1.0, 1, 0, 2.0, nlower, norbs, 0, pp, 5);
    if (!(near(pp[0], 2.0) && near(pp[1], 0.0) && near(pp[2], 0.0))) {
        std::printf("FAIL syrk: %g %g %g\n", pp[0], pp[1], pp[2]);
        ok = false;
    }

    // --- Non-identity: c = [[1,0],[0,1]] scaled columns; ndubl=2, nsingl=0 ---
    // P = 2*c*c^T = 2*I -> packed (2,0,2)
    double c2[4] = {0.6, 0.8, -0.8, 0.6};  // orthogonal columns
    for (int i = 0; i < 3; ++i) pp[i] = -9.0;
    density_for_GPU(c2, 1.0, 2, 0, 2.0, nlower, norbs, 0, pp, 3);
    // col1=(0.6,0.8), col2=(-0.8,0.6): c*c^T = I (orthonormal)
    if (!(near(pp[0], 2.0) && near(pp[1], 0.0) && near(pp[2], 2.0))) {
        std::printf("FAIL orthonormal dgemm: %g %g %g\n", pp[0], pp[1], pp[2]);
        ok = false;
    }

    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
