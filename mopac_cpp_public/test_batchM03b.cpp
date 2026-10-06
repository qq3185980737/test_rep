// test_batchM03b.cpp — M03 两电子/重叠积分核心 语义断言
// 覆盖：bfn(B 积分三分支)、quadr(二次插值)、epseta(机器常量)、xxx(标签构造)
#include <cstdio>
#include <cmath>
#include <string>
#include "bfn.h"
#include "quadr.h"
#include "epseta.h"
#include "xxx.h"
#include "overlaps_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_bfn() {
    std::printf("Test 1: bfn B integrals (three branches)\n");
    // set factorial table fact[0..17] used by series branch
    overlaps_C::fact[0] = 1.0;
    for (int i = 1; i <= 17; ++i) overlaps_C::fact[i] = overlaps_C::fact[i-1] * i;
    // reference values computed independently from the F90 formula
    double bf[14];
    bfn(0.0001, bf);   // small-x branch (label 90)
    CHK_D(bf[1], 2.000000003, 1e-7, "bfn small-x bf1");
    CHK_D(bf[2], -6.666666673e-05, 1e-9, "bfn small-x bf2");
    CHK_D(bf[3], 0.6666666687, 1e-7, "bfn small-x bf3");
    bfn(1.5, bf);      // series branch (last=12)
    CHK_D(bf[1], 2.839039273, 1e-7, "bfn series bf1");
    CHK_D(bf[2], -1.2438533, 1e-7, "bfn series bf2");
    CHK_D(bf[4], -0.7754097489, 1e-7, "bfn series bf4");
    bfn(3.5, bf);      // recurrence branch (|x|>3)
    CHK_D(bf[1], 9.452929879, 1e-7, "bfn recurrence bf1");
    CHK_D(bf[2], -6.769348418, 1e-7, "bfn recurrence bf2");
    CHK_D(bf[5], 4.100617391, 1e-7, "bfn recurrence bf5");
}

static void test_quadr() {
    std::printf("Test 2: quadr quadratic through 3 points\n");
    double a = 0, b = 0, c = 0;
    // f0=1 at 0, f1=4 at 1, f2=9 at 2 => f = 1 + 2x + x^2
    quadr(1.0, 4.0, 9.0, 1.0, 2.0, a, b, c);
    CHK_D(a, 1.0, 1e-12, "quadr a=f0");
    CHK_D(b, 2.0, 1e-12, "quadr b");
    CHK_D(c, 1.0, 1e-12, "quadr c");
    // asymmetric: f0=0 at 0, f1=1 at 1, f2=8 at 4  =>  f = x^2/2 + x/2?
    // solve: c*1+b*1=1; c*16+b*4=8 -> 15c+3b=7... use formula directly:
    quadr(0.0, 1.0, 8.0, 1.0, 4.0, a, b, c);
    // c = (4*(1-0) - 1*(8-0))/(4*1*(1-4)) = (4-8)/(-12) = 1/3
    // b = (1-0 - (1/3)*1)/1 = 2/3
    CHK_D(a, 0.0, 1e-12, "quadr2 a");
    CHK_D(b, 2.0/3.0, 1e-12, "quadr2 b");
    CHK_D(c, 1.0/3.0, 1e-12, "quadr2 c");
}

static void test_epseta() {
    std::printf("Test 3: epseta machine constants\n");
    double eps = 0, eta = 0;
    epseta(eps, eta);
    // eps = max(Tiny, 1e-39) = 1e-39; eta = max(Epsilon, 1e-17) = DBL_EPSILON
    CHK_D(eps, 1e-39, 1e-50, "epseta eps = 1e-39");
    CHK_D(eta, 2.220446049250313e-16, 1e-30, "epseta eta = DBL_EPSILON");
}

static void test_xxx() {
    std::printf("Test 4: xxx connectivity label\n");
    std::string r;
    xxx('R', 1, 2, 0, 3, r);
    CHECK(r == "R123", "xxx label R123");
    xxx('P', 12, 3, 0, 0, r);
    CHECK(r == "P123", "xxx label P123 (12 -> '12')");
    xxx('F', 0, 0, 0, 0, r);
    CHECK(r == "F", "xxx label F only");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M03b batch — two-electron / overlap integral core\n");
    test_bfn();
    test_quadr();
    test_epseta();
    test_xxx();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
