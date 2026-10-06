// test_batchM06a.cpp — M06 梯度：deri0(超级矩阵对角+标量) 数值 + deri21(PCA基) 结构
#include <cstdio>
#include <cmath>
#include <vector>
#include "deri0.h"
#include "deri21.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_deri0() {
    std::printf("Test 1: deri0 super-matrix diagonal + scalar\n");
    // n=6, closed 2 / open 2 / virtual 2, fract=0.5, e = {-10,-8,-6,-4,-2,0}
    std::vector<double> e(7, 0.0);
    for (int i = 1; i <= 6; ++i) e[i] = -12.0 + 2.0 * i;  // -10,-8,-6,-4,-2,0
    int nbo_arr[4] = {0, 2, 2, 2};
    std::vector<int> nbo(nbo_arr, nbo_arr + 4);
    std::vector<double> diag(14, -1.0), scalar(14, -1.0);
    deri0(e, 6, scalar, diag, 0.5, nbo);
    // hand-computed from F90 formula:
    // open-closed: (e3-e1)/1.501=2.664890, (e4-e1)/1.501=3.997335,
    //              (e3-e2)/1.501=1.332445, (e4-e2)/1.501=2.664890
    CHK_D(diag[1], 4.0 / 1.501, 1e-6, "deri0 oc1");
    CHK_D(diag[2], 6.0 / 1.501, 1e-6, "deri0 oc2");
    CHK_D(diag[3], 2.0 / 1.501, 1e-6, "deri0 oc3");
    // virtual-closed: (e5-e1)/2=4, (e6-e1)/2=5, (e5-e2)/2=3, (e6-e2)/2=4
    CHK_D(diag[5], 4.0, 1e-12, "deri0 vc1");
    CHK_D(diag[6], 5.0, 1e-12, "deri0 vc2");
    // virtual-open: (e5-e3)/0.501=7.984032, (e6-e4)/0.501=7.984032
    CHK_D(diag[9], 4.0 / 0.501, 1e-6, "deri0 vo1");
    CHK_D(diag[12], 4.0 / 0.501, 1e-6, "deri0 vo4");
    // scalar[i] = sqrt(1/max(0.3*diag, diag-2.36))
    CHK_D(scalar[1], std::sqrt(1.0 / std::max(0.3 * (4.0 / 1.501), (4.0 / 1.501) - 2.36)), 1e-9, "deri0 scalar1");
    CHK_D(scalar[5], std::sqrt(1.0 / std::max(1.2, 1.64)), 1e-9, "deri0 scalar5");
    CHK_D(scalar[12], std::sqrt(1.0 / std::max(0.3 * (4.0 / 0.501), (4.0 / 0.501) - 2.36)), 1e-9, "deri0 scalar12");
}

static void test_deri21() {
    std::printf("Test 2: deri21 PCA orthonormal basis\n");
    // a: 4x3 column-major: col1=(1,0,0,1), col2=(0,1,0,1), col3=(0,0,1,1)
    double a[12] = {1,0,0,1,  0,1,0,1,  0,0,1,1};
    double vnert[9] = {}, pnert[3] = {};
    double b[12] = {};
    int ncut = 0;
    deri21(a, 3, 4, 2, vnert, pnert, b, ncut);
    CHECK(ncut == 2, "deri21 ncut==2");
    CHK_D(pnert[0], 2.0, 1e-6, "deri21 pnert1=sqrt(4)");
    CHK_D(pnert[1], 1.0, 1e-6, "deri21 pnert2=sqrt(1)");
    // B col1 = A*(1,1,1)/(2*sqrt3) = (1,1,1,3)/(2sqrt3)
    double s23 = 1.0 / (2.0 * std::sqrt(3.0));
    CHK_D(b[0], s23, 1e-6, "deri21 b11");
    CHK_D(b[1], s23, 1e-6, "deri21 b21");
    CHK_D(b[3], 3.0 * s23, 1e-6, "deri21 b41");
    // orthonormality: B'B = I2
    double b11 = 0, b12 = 0, b22 = 0;
    for (int i = 0; i < 4; ++i) {
        b11 += b[i] * b[i];
        b12 += b[i] * b[4 + i];
        b22 += b[4 + i] * b[4 + i];
    }
    CHK_D(b11, 1.0, 1e-6, "deri21 b col1 norm");
    CHK_D(b22, 1.0, 1e-6, "deri21 b col2 norm");
    CHK_D(b12, 0.0, 1e-6, "deri21 b cols orthogonal");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06a batch — deri0 / deri21\n");
    test_deri0();
    test_deri21();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
