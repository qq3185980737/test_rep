// test_batchM03b2.cpp — M03 两电子族：fock2 分支数值 + dhc/rotate 结构
// 注：rotatd/elenuc/nddo_to_point/diat 的 F90 源码不在 MOPAC 2016 树
//（老 F77 NDDO 积分内核遗留），此处提供测试 stub 仅用于验证 C++ 结构
// 编译与 Fock 索引逻辑；积分内核数值对拍待外部内核接入。
#include <cstdio>
#include <cmath>
#include <vector>
#include "fock2.h"
#include "dhc.h"
#include "rotate.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "parameters_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

extern "C" void rotatd_(int*, int*, const double*, const double*, double* w, int* kr, double* enuc) {
    *kr = 1; *enuc = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
    for (int i = 0; i < 171; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double*, double*, double*, double*, double*, int*, int*) {}
extern "C" void diat_(int*, int*, double*, double* smat) {
    for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}

static void test_fock2_light_light() {
    std::printf("Test 1: fock2 light-light (1-orbital atoms)\n");
    molkst_C::numcal = 1;
    molkst_C::id = 0;
    molkst_C::norbs = 2;
    std::vector<int> nfirst(3, 0), nlast(3, 0);
    nfirst[1] = 1; nlast[1] = 1;
    nfirst[2] = 2; nlast[2] = 2;
    // packed 1-based: p1=P11, p2=P12, p3=P22
    std::vector<double> f(4, 0.0), ptot(4, 0.0), p(4, 0.0);
    ptot[1] = 2.0; ptot[2] = 0.5; ptot[3] = 3.0;
    p[1] = 2.0; p[2] = 0.5; p[3] = 3.0;
    std::vector<double> w(2, 0.0), wj(2, 0.0), wk(2, 0.0);
    w[1] = 0.3;
    fock2(f, ptot, p, w, wj, wk, 2, nfirst, nlast, 1);
    // i1=i1fact(1)=1, j1=i1fact(2)=3, ij=1+2-1=2
    // f1 += ptot3*w1 = 0.9 ; f3 += ptot1*w1 = 0.6 ; f2 -= p2*w1 = -0.15
    CHK_D(f[1], 0.9, 1e-12, "fock2 l-l f1");
    CHK_D(f[2], -0.15, 1e-12, "fock2 l-l f2");
    CHK_D(f[3], 0.6, 1e-12, "fock2 l-l f3");
}

static void test_fock2_onecenter() {
    std::printf("Test 2: fock2 mode=2 fock1dorbs (single atom, 2 orbitals)\n");
    molkst_C::numcal = 2;
    molkst_C::id = 0;
    molkst_C::norbs = 2;
    std::vector<int> nfirst(2, 0), nlast(2, 0);
    nfirst[1] = 1; nlast[1] = 2;
    std::vector<double> f(4, 0.0), ptot(4, 0.0), p(4, 0.0);
    ptot[1] = 2.0; ptot[2] = 0.0; ptot[3] = 5.0;
    p = ptot; p[1] = 1.0; p[3] = 4.0;
    std::vector<double> w(10, 0.0), wj(10, 0.0), wk(10, 0.0);
    for (int i = 1; i <= 9; ++i) w[i] = 1.0;
    fock2(f, ptot, p, w, wj, wk, 1, nfirst, nlast, 2);
    // wloc all ones => f[ij] = sum_k,l (ptot[ijp]-pa[ijp]) per F90 fock1dorbs
    // f1 = (2-1)+(5-4)=2 ; f2 same by symmetry of contributions =2 ; f3 = 2
    CHK_D(f[1], 2.0, 1e-12, "fock1dorbs f1");
    CHK_D(f[2], 2.0, 1e-12, "fock1dorbs f2");
    CHK_D(f[3], 2.0, 1e-12, "fock1dorbs f3");
}

static void test_rotate_layout() {
    std::printf("Test 3: rotate 1-based w layout (stub kernel)\n");
    molkst_C::numcal = 3;
    molkst_C::id = 0;
    parameters_C::natorb[1] = 1; parameters_C::natorb[6] = 1;
    double w[2027] = {};
    for (int i = 1; i < 2026; ++i) w[i] = -1.0;  // sentinel: untouched positions
    double e1b[45] = {}, e2a[45] = {};
    int kr = 0;
    double enuc = 0.0;
    double xi[3] = {0, 0, 0}, xj[3] = {2.0, 0, 0};
    rotate(6, 1, xi, xj, w, kr, e1b, e2a, enuc);
    CHECK(kr == 1, "rotate kr==1 (stub kernel)");
    CHK_D(enuc, 0.0, 1e-15, "rotate enuc==0");
    CHK_D(w[1], -1.0, 1e-15, "rotate w[1] untouched (1-based)");
    CHK_D(w[2025], -1.0, 1e-15, "rotate w[2025] untouched");
}

static void test_dhc_molecule() {
    std::printf("Test 4: dhc structural run (molecular id=0, stub kernels)\n");
    molkst_C::numcal = 4;
    molkst_C::id = 0;
    molkst_C::uhf = false;
    molkst_C::cutofp = 100.0;
    MOZYME_C::cutofs = 100.0;
    molkst_C::norbs = 2;
    double p[10] = {}, pa[10] = {}, pb[10] = {};
    pa[0] = 0.5; pb[0] = 0.5; p[0] = 1.0;
    // two atoms 2 A apart, 1 orbital each
    double xi[6] = {0, 0, 0, 2.0, 0, 0};
    int nat[2] = {1, 1};
    double dener = 99.0;
    dhc(p, pa, pb, xi, nat, 1, 1, 1, 1, dener, 1);
    // h1elec: diat stub -> shmat=0 -> h=f=0; rotate: kernel stub -> w=0, e1b/e2a=0
    // fock2 l-l on zero w -> f=0; helect on zero h/f -> 0; dener = 0
    CHK_D(dener, 0.0, 1e-12, "dhc dener==0 with stub kernels");
    // non-first-call path (icalcn==numcal) must not disturb state
    dhc(p, pa, pb, xi, nat, 1, 1, 1, 1, dener, 1);
    CHK_D(dener, 0.0, 1e-12, "dhc second call stable");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M03b2 batch — fock2 branches + dhc/rotate structure\n");
    test_fock2_light_light();
    test_fock2_onecenter();
    test_rotate_layout();
    test_dhc_molecule();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
