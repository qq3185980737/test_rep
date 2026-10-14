// test_batchM08a.cpp  M08 eigenvectors_LAPACK (dsyevd-equivalent) numerics
#include <cstdio>
#include <cmath>
#include <vector>
#include "eigenvectors_LAPACK.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_2x2() {
    std::printf("Test 1: 2x2 [[2,1],[1,2]] -> eigvals 1,3\n");
    double xmat[3] = {2.0, 1.0, 2.0};  // packed upper, 0-based
    double vecs[4] = {}, vals[2] = {};
    eigenvectors_LAPACK(vecs, xmat, vals, 2);
    CHK_D(vals[0], 1.0, 1e-9, "eigval0 1");
    CHK_D(vals[1], 3.0, 1e-9, "eigval1 3");
    // column-major: col0 = (1,-1)/sqrt2, col1 = (1,1)/sqrt2
    CHK_D(std::fabs(vecs[0]) + std::fabs(vecs[1]), std::sqrt(2.0), 1e-9,
          "col0 norm");
    // orthogonality
    double dot = vecs[0] * vecs[2] + vecs[1] * vecs[3];
    CHK_D(dot, 0.0, 1e-9, "cols orthogonal");
}

static void test_3x3_diag_order() {
    std::printf("Test 2: 3x3 diag [3,1,2] -> sorted ascending 1,2,3\n");
    double xmat[6] = {3.0, 0.0, 1.0, 0.0, 0.0, 2.0};
    double vecs[9] = {}, vals[3] = {};
    eigenvectors_LAPACK(vecs, xmat, vals, 3);
    CHK_D(vals[0], 1.0, 1e-6, "sorted eigval0 1");
    CHK_D(vals[1], 2.0, 1e-6, "sorted eigval1 2");
    CHK_D(vals[2], 3.0, 1e-6, "sorted eigval2 3");
    // eigenvector of eigval 1 (idx 1): column = e1 (in col-major: vecs[1])
    CHK_D(vecs[1], 1.0, 1e-6, "eigvec col0 = e1 (sorted)");
}

static void test_3x3_general() {
    std::printf("Test 3: 3x3 [[2,1,0],[1,3,1],[0,1,2]] -> 1,2,4\n");
    double xmat[6] = {2.0, 1.0, 3.0, 0.0, 1.0, 2.0};
    // packed: (0,0)=2,(0,1)=1,(1,1)=3,(0,2)=0,(1,2)=1,(2,2)=2
    double vecs[9] = {}, vals[3] = {};
    eigenvectors_LAPACK(vecs, xmat, vals, 3);
    CHK_D(vals[0], 1.0, 1e-8, "eigval0 1");
    CHK_D(vals[1], 2.0, 1e-8, "eigval1 2");
    CHK_D(vals[2], 4.0, 1e-8, "eigval2 4");
    // eigenvector columns orthonormal
    for (int c = 0; c < 3; ++c) {
        double n2 = 0.0;
        for (int r = 0; r < 3; ++r) n2 += vecs[r + c * 3] * vecs[r + c * 3];
        CHK_D(n2, 1.0, 1e-8, "col orthonormal");
    }
    // residual: A*v - lam*v ~ 0 (rebuild A)
    double A[9] = {2,1,0, 1,3,1, 0,1,2};
    for (int c = 0; c < 3; ++c) {
        double res = 0.0;
        for (int r = 0; r < 3; ++r) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += A[r * 3 + k] * vecs[k + c * 3];
            s -= vals[c] * vecs[r + c * 3];
            res += s * s;
        }
        CHK_D(res, 0.0, 1e-6, "residual col");
    }
}

static void test_degenerate_split() {
    std::printf("Test 4: degeneracy split by i*1e-10 perturbation\n");
    double xmat[3] = {1.0, 0.0, 1.0};  // identity 2x2
    double vecs[4] = {}, vals[2] = {};
    eigenvectors_LAPACK(vecs, xmat, vals, 2);
    // split 1+1e-10 vs 1+2e-10 -> strictly increasing
    CHECK(vals[0] < vals[1], "degenerate split ascending");
    CHK_D(vals[0], 1.0 + 1.0e-10, 1e-12, "perturbation level 1e-10");
    CHK_D(vals[1], 1.0 + 2.0e-10, 1e-12, "perturbation level 2e-10");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08a batch - eigenvectors_LAPACK\n");
    test_2x2();
    test_3x3_diag_order();
    test_3x3_general();
    test_degenerate_split();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
