// test_batchM06b.cpp — M06：dfock2 两电子 Fock 导数（轻-轻/轻-重数值 + 重-重结构）
#include <cstdio>
#include <cmath>
#include <vector>
#include "dfock2.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_dfock2_light_light() {
    std::printf("Test 1: dfock2 light-light\n");
    molkst_C::numcal = 1;
    molkst_C::norbs = 2;
    int nfirst[3] = {0, 1, 2}, nlast[3] = {0, 1, 2};
    double f[10] = {}, ptot[10] = {}, p[10] = {}, w[10] = {};
    ptot[1] = 2.0; ptot[2] = 0.5; ptot[3] = 3.0;
    p[1] = 2.0; p[2] = 0.5; p[3] = 3.0;
    w[1] = 0.4;
    dfock2(f, ptot, p, w, 2, nfirst, nlast, 1);
    // i1=i1fact(1)=1, j1=i1fact(2)=3, ij=3+1-2=2
    // f1 += ptot3*w1 = 1.2 ; f3 += ptot1*w1 = 0.8 ; f2 -= p2*w1 = -0.2
    CHK_D(f[1], 1.2, 1e-12, "dfock2 l-l f1");
    CHK_D(f[2], -0.2, 1e-12, "dfock2 l-l f2");
    CHK_D(f[3], 0.8, 1e-12, "dfock2 l-l f3");
}

static void test_dfock2_heavy_light() {
    std::printf("Test 2: dfock2 light atom ii (1 orb) - heavy atom jj (4 orbs)\n");
    molkst_C::numcal = 2;
    molkst_C::norbs = 5;
    int nfirst[3] = {0, 1, 2}, nlast[3] = {0, 1, 5};
    double f[30] = {}, ptot[30] = {}, p[30] = {}, w[30] = {};
    ptot[1] = 1.0;  // P11 (ll = i1fact(1) = 1)
    for (int i = 1; i <= 10; ++i) w[i] = 0.4;
    dfock2(f, ptot, p, w, 2, nfirst, nlast, 1);
    // Coulomb: heavy-atom diagonal j1 sequence (ja=2): i=0->3, i=1->6, i=2->10, i=3->15
    CHK_D(f[3], 0.4, 1e-12, "dfock2 h-l f3");
    CHK_D(f[6], 0.4, 1e-12, "dfock2 h-l f6");
    CHK_D(f[10], 0.4, 1e-12, "dfock2 h-l f10");
    CHK_D(f[15], 0.4, 1e-12, "dfock2 h-l f15");
    CHK_D(f[5], 0.4, 1e-12, "dfock2 h-l f5 offdiag");
    CHK_D(f[7], 0.0, 1e-12, "dfock2 h-l f7 unwritten");
    // sumoff/sumdia: ptot[j1] all zero => f[1] unchanged (0)
    CHK_D(f[1], 0.0, 1e-12, "dfock2 h-l f1 (sum zero)");
    // exchange: p all zero => no change
    CHK_D(f[2], 0.0, 1e-12, "dfock2 h-l f2 (exchange zero)");
}

static void test_dfock2_heavy_heavy() {
    std::printf("Test 3: dfock2 heavy-heavy structural (jab+kab)\n");
    molkst_C::numcal = 3;
    molkst_C::norbs = 8;
    int nfirst[3] = {0, 1, 5}, nlast[3] = {0, 4, 8};
    double f[80] = {}, ptot[80] = {}, p[80] = {}, w[120] = {};
    for (int i = 1; i <= 3; ++i) ptot[i] = 1.0;  // atom-1 block density
    for (int i = 1; i <= 100; ++i) w[i] = 0.1;
    dfock2(f, ptot, p, w, 2, nfirst, nlast, 1);
    bool changed = false;
    for (int i = 1; i <= 36; ++i) if (std::fabs(f[i]) > 1e-12) changed = true;
    CHECK(changed, "dfock2 h-h f updated (jab/kab)");
}

static void test_dfock2_case2_like() {
    std::printf("Test 4: dfock2 second-call stability (itype save path)\n");
    molkst_C::numcal = 4;
    molkst_C::norbs = 2;
    int nfirst[3] = {0, 1, 2}, nlast[3] = {0, 1, 2};
    double f[10] = {}, ptot[10] = {}, p[10] = {}, w[10] = {};
    ptot[1] = 1.0; ptot[2] = 0.0; ptot[3] = 1.0;
    p[1] = 1.0; p[2] = 0.0; p[3] = 1.0;
    w[1] = 0.5;
    dfock2(f, ptot, p, w, 2, nfirst, nlast, 1);
    CHK_D(f[1], 0.5, 1e-12, "dfock2 2nd-call f1 (fresh icalcn)");
    CHK_D(f[3], 0.5, 1e-12, "dfock2 2nd-call f3");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06b batch — dfock2\n");
    test_dfock2_light_light();
    test_dfock2_heavy_light();
    test_dfock2_heavy_heavy();
    test_dfock2_case2_like();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
