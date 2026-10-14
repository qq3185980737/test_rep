// test_batchM12.cpp — M12 线性代数/数值工具 语义断言（数值对拍）
// 覆盖：minv/osinv(求逆 A*A^-1≈I)、mxm(A*B)、mxmt(A*B^T)、mtxm(A^T*B)、
//       mxv(A*x)、mtxmc(下三角打包)、mult(C*S)、dot、sort、supdot、mat33、rsp、linpack
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include "minv.h"
#include "osinv.h"
#include "mxm.h"
#include "mxmt.h"
#include "mtxm.h"
#include "mxv.h"
#include "mtxmc.h"
#include "mult.h"
#include "dot.h"
#include "sort.h"
#include "supdot.h"
#include "mat33.h"
#include "rsp.h"
#include "linpack.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)

// 3x3 column-major mul
static void mul33(double* r, const double* a, const double* b) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0;
            for (int k = 0; k < 3; ++k) s += a[k * 3 + i] * b[j * 3 + k];
            r[j * 3 + i] = s;
        }
}

static void test_minv() {
    std::printf("Test 1: minv Gauss-Jordan inverse\n");
    // A = [4 7 2; 3 6 1; 2 3 5] column-major
    double a[9] = {4, 3, 2, 7, 6, 3, 2, 1, 5};
    double a0[9];
    for (int i = 0; i < 9; ++i) a0[i] = a[i];
    double d = 0;
    minv(a, 3, d);
    // det = 4*(30-3) - 7*(15-2) + 2*(9-12) = 108 - 91 - 6 = 11
    CHECK(std::fabs(d - 11.0) < 1e-9, "minv determinant = 11");
    double prod[9];
    mul33(prod, a0, a);  // A * A^-1
    bool ok = true;
    for (int i = 0; i < 9; ++i) {
        double want = (i / 3 == i % 3) ? 1.0 : 0.0;
        if (std::fabs(prod[i] - want) > 1e-9) ok = false;
    }
    CHECK(ok, "minv A * A^-1 = I");
}

static void test_osinv() {
    std::printf("Test 2: osinv symmetric inverse\n");
    double a[9] = {4, 3, 2, 3, 6, 3, 2, 3, 5};  // symmetric, det=4*(30-9)-3*(15-6)+2*(9-12)=84-27-6=51
    double a0[9];
    for (int i = 0; i < 9; ++i) a0[i] = a[i];
    double d = 0;
    osinv(a, 3, d);
    double prod[9];
    mul33(prod, a0, a);
    bool ok = true;
    for (int i = 0; i < 9; ++i) {
        double want = (i / 3 == i % 3) ? 1.0 : 0.0;
        if (std::fabs(prod[i] - want) > 1e-9) ok = false;
    }
    CHECK(ok, "osinv A * A^-1 = I");
}

static void test_mxm() {
    std::printf("Test 3: mxm A*B\n");
    // A=[1 2;3 4] col-major {1,3,2,4}, B=[5 6;7 8] col-major {5,7,6,8}
    double a[4] = {1, 3, 2, 4}, b[4] = {5, 7, 6, 8}, c[4];
    mxm(a, 2, b, 2, c, 2);
    // A*B = [19 22; 43 50]
    CHECK(std::fabs(c[0] - 19) < 1e-9 && std::fabs(c[1] - 43) < 1e-9 &&
          std::fabs(c[2] - 22) < 1e-9 && std::fabs(c[3] - 50) < 1e-9,
          "mxm 2x2 product");
    // rectangular 2x3 * 3x2
    double a2[6] = {1, 4, 2, 5, 3, 6};           // [1 2 3; 4 5 6]
    double b2[6] = {7, 9, 11, 8, 10, 12};        // [7 8; 9 10; 11 12] col-major 3x2
    double c2[4];
    mxm(a2, 2, b2, 3, c2, 2);
    // = [1*7+2*9+3*11, 1*8+2*10+3*12; 4*7+5*9+6*11, 4*8+5*10+6*12]
    // = [58, 64; 139, 154]
    CHECK(std::fabs(c2[0] - 58) < 1e-9 && std::fabs(c2[1] - 139) < 1e-9 &&
          std::fabs(c2[2] - 64) < 1e-9 && std::fabs(c2[3] - 154) < 1e-9,
          "mxm rectangular product");
}

static void test_mxmt() {
    std::printf("Test 4: mxmt A*B^T\n");
    // A=[1 2;3 4] col {1,3,2,4}; B = [5 6;7 8] col {5,7,6,8}; B^T = [5 7;6 8]
    double a[4] = {1, 3, 2, 4}, b[4] = {5, 7, 6, 8}, c[4];
    mxmt(a, 2, b, 2, c, 2);
    // A*B^T = [1*5+2*6, 1*7+2*8; 3*5+4*6, 3*7+4*8] = [17,23;39,53]
    CHECK(std::fabs(c[0] - 17) < 1e-9 && std::fabs(c[1] - 39) < 1e-9 &&
          std::fabs(c[2] - 23) < 1e-9 && std::fabs(c[3] - 53) < 1e-9,
          "mxmt A*B^T");
}

static void test_mtxm() {
    std::printf("Test 5: mtxm A^T*B\n");
    // A = [1 2;3 4] col {1,3,2,4} (A^T = [1 3;2 4])
    // B = [5 6;7 8] col {5,7,6,8}
    double a[4] = {1, 3, 2, 4}, b[4] = {5, 7, 6, 8}, c[4];
    mtxm(a, 2, b, 2, c, 2);
    // A^T*B = [1*5+3*7, 1*6+3*8; 2*5+4*7, 2*6+4*8] = [26,30;38,44]
    CHECK(std::fabs(c[0] - 26) < 1e-9 && std::fabs(c[1] - 38) < 1e-9 &&
          std::fabs(c[2] - 30) < 1e-9 && std::fabs(c[3] - 44) < 1e-9,
          "mtxm A^T*B");
}

static void test_mxv() {
    std::printf("Test 6: mxv A*x (vectors 1-based)\n");
    double a[4] = {1, 3, 2, 4};  // [1 2;3 4]
    // vecx 1-based: index 1..2
    std::vector<double> vx(3, 0.0), vy(3, 0.0);
    vx[1] = 5; vx[2] = 7;
    mxv(a, 2, &vx[0], 2, &vy[0]);
    // A*x = [1*5+2*7; 3*5+4*7] = [19;43]
    CHECK(std::fabs(vy[1] - 19) < 1e-9 && std::fabs(vy[2] - 43) < 1e-9,
          "mxv A*x");
}

static void test_mtxmc() {
    std::printf("Test 7: mtxmc lower triangle\n");
    // A = B = [1 2;3 4] col {1,3,2,4} (nbr=2 rows, nar=2 cols)
    double a[4] = {1, 3, 2, 4}, b[4] = {1, 3, 2, 4}, c[3];
    mtxmc(a, 2, b, 2, c);
    // C = A^T*B = [1*1+3*3, 1*2+3*4; 2*1+4*3, 2*2+4*4] = [10,14;14,20]
    // packed lower triangle canonical: c = {c(1,1), c(2,1), c(2,2)} = {10,14,20}
    CHECK(std::fabs(c[0] - 10) < 1e-9 && std::fabs(c[1] - 14) < 1e-9 &&
          std::fabs(c[2] - 20) < 1e-9, "mtxmc packed lower triangle");
}

static void test_mult() {
    std::printf("Test 8: mult C*S -> S*C (Mulliken back-transform)\n");
    double c[4] = {1, 3, 2, 4};  // [1 2;3 4]
    double s[4] = {5, 7, 6, 8};  // [5 6;7 8]
    double v[4];
    mult(c, s, v, 2);
    // F90: vecs(j,i) = sum_k c(k,i)*s(j,k) => vecs = S*C
    // S*C = [5*1+6*3, 5*2+6*4; 7*1+8*3, 7*2+8*4] = [23,31;34,46] col {23,31,34,46}
    CHECK(std::fabs(v[0] - 23) < 1e-9 && std::fabs(v[1] - 31) < 1e-9 &&
          std::fabs(v[2] - 34) < 1e-9 && std::fabs(v[3] - 46) < 1e-9,
          "mult S*C");
}

static void test_dot() {
    std::printf("Test 9: dot product (1-based vectors)\n");
    std::vector<double> x(3, 0.0), y(3, 0.0);
    x[1] = 1.5; x[2] = -2.0;
    y[1] = 4.0; y[2] = 0.5;
    CHECK(std::fabs(dot(x, y, 2) - (6.0 - 1.0)) < 1e-12, "dot = 5.0");
}

static void test_sort() {
    std::printf("Test 10: sort selection (real vals + complex vec)\n");
    float val[5] = {0, 3.0f, 1.0f, 4.0f, 2.0f};   // 1-based values
    std::complex<float> vec[16];                   // 4x4 col-major (n rows x n cols), 0-based
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            vec[j * 4 + i] = std::complex<float>((float)(j + 1), (float)i);
    // val(1..4) = {3,1,4,2}; after sort -> {1,2,3,4} ascending by selection
    sort(val, vec, 4);
    bool ok = std::fabs(val[1] - 1.0f) < 1e-5 && std::fabs(val[2] - 2.0f) < 1e-5 &&
              std::fabs(val[3] - 3.0f) < 1e-5 && std::fabs(val[4] - 4.0f) < 1e-5;
    CHECK(ok, "sort values ascending");
    // val column j pairs swap with vec columns: val(1..4)={3,1,4,2}
    // sel sort: i=1 swaps col2<->col1, i=2 swaps col4<->col2, i=3 swaps col4<->col3
    // col1 ends up with original col2 = (2,0),(2,1),(2,2),(2,3)
    // col2 ends up with original col4 = (4,0)...
    CHECK(std::fabs(vec[0].real() - 2.0f) < 1e-5 &&
          std::fabs(vec[4].real() - 4.0f) < 1e-5, "sort vec col swap");
}

static void test_supdot() {
    std::printf("Test 11: supdot S = H*G (H symmetric, packed lower)\n");
    // H symmetric [4 2 1; 2 5 3; 1 3 6], packed lower canonical:
    // h[1]=H11=4, h[2]=H21=2, h[3]=H22=5, h[4]=H31=1, h[5]=H32=3, h[6]=H33=6
    // G = (1, -2, 3); S = H*G = (4*1+2*(-2)+1*3, 2*1+5*(-2)+3*3, 1*1+3*(-2)+6*3)
    //   = (4-4+3, 2-10+9, 1-6+18) = (3, 1, 13)
    std::vector<double> h(7, 0.0), g(4, 0.0), s(4, 0.0);
    h[1] = 4; h[2] = 2; h[3] = 5; h[4] = 1; h[5] = 3; h[6] = 6;
    g[1] = 1; g[2] = -2; g[3] = 3;
    supdot(&s[0], &h[0], &g[0], 3);
    CHECK(std::fabs(s[1] - 3) < 1e-12 && std::fabs(s[2] - 1) < 1e-12 &&
          std::fabs(s[3] - 13) < 1e-12, "supdot H*G");
}

static void test_mat33() {
    std::printf("Test 12: mat33 A^T B A (3x3, 1-based a[1..9])\n");
    // A = diag(1,2,3) col-major 1-based: a[1]=1,a[5]=2,a[9]=3
    // B = all-ones
    double a[10] = {0, 1, 0, 0, 0, 2, 0, 0, 0, 3};
    double b[10] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    double c[10] = {0};
    mat33(a, b, c);
    // t = B*A = each col sum of A -> rows [1,2,3]
    // C = A^T*t = diag(1,2,3)*t -> rows [1,2,3],[2,4,6],[3,6,9]
    // col-major 1-based: c[1..9] = {1,2,3, 2,4,6, 3,6,9}
    double want[10] = {0, 1, 2, 3, 2, 4, 6, 3, 6, 9};
    bool ok = true;
    for (int i = 1; i <= 9; ++i)
        if (std::fabs(c[i] - want[i]) > 1e-9) ok = false;
    CHECK(ok, "mat33 A^T B A");
}

static void test_rsp() {
    std::printf("Test 13: rsp symmetric eigenvalues\n");
    // A = [2 1 0; 1 2 1; 0 1 2] eigenvalues 2+sqrt(2), 2, 2-sqrt(2)
    // packed lower 0-based: a[0]=A11=2, a[1]=A21=1, a[2]=A22=2, a[3]=A31=0, a[4]=A32=1, a[5]=A33=2
    double a[6] = {2, 1, 2, 0, 1, 2};
    double root[4], vect[12];
    rsp(a, 3, root, vect);
    double e0 = 2.0 + std::sqrt(2.0), e1 = 2.0, e2 = 2.0 - std::sqrt(2.0);
    // ascending: e2, e1, e0
    CHECK(std::fabs(root[0] - e2) < 1e-9 && std::fabs(root[1] - e1) < 1e-9 &&
          std::fabs(root[2] - e0) < 1e-9, "rsp eigenvalues ascending");
}

static void test_linpack() {
    std::printf("Test 14: linpack dgefa/dgedi\n");
    // A = [2 1 1; 4 -6 0; -2 7 2] det = 2*(-12) - 1*8 + 1*(28-12) = -24-8+16=-16
    double a[9] = {2, 4, -2, 1, -6, 7, 1, 0, 2};
    double a0[9];
    for (int i = 0; i < 9; ++i) a0[i] = a[i];
    int ipvt[4], info = 0;   // 1-based logical: ipvt[1..n] -> need n+1
    dgefa(a, 3, 3, ipvt, info);
    CHECK(info == 0, "dgefa no singular");
    double det[2], work[4];  // work[i] i=1..n
    dgedi(a, 3, 3, ipvt, det, work, 10);   // job=10: det only
    // det represented as det[0]*10^det[1]; det = -16 -> det[0]=-1.6, det[1]=1
    CHECK(std::fabs(det[0] * std::pow(10.0, det[1]) - (-16.0)) < 1e-9,
          "dgedi determinant");
    // job=1: inverse only
    for (int i = 0; i < 9; ++i) a[i] = a0[i];
    dgefa(a, 3, 3, ipvt, info);
    dgedi(a, 3, 3, ipvt, det, work, 1);
    double prod[9];
    mul33(prod, a0, a);
    bool ok = true;
    for (int i = 0; i < 9; ++i) {
        double want = (i / 3 == i % 3) ? 1.0 : 0.0;
        if (std::fabs(prod[i] - want) > 1e-9) ok = false;
    }
    CHECK(ok, "dgedi inverse A*A^-1=I");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M12 batch — linear algebra / numeric tools\n");
    test_minv();
    test_osinv();
    test_mxm();
    test_mxmt();
    test_mtxm();
    test_mxv();
    test_mtxmc();
    test_mult();
    test_dot();
    test_sort();
    test_supdot();
    test_mat33();
    test_rsp();
    test_linpack();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
