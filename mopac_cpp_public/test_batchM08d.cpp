// test_batchM08d.cpp  M08 dijkl1 (atom-fixed 2e integrals) + dijkl2 (relaxation)
#include <cstdio>
#include <cmath>
#include <vector>
#include "dijkl1.h"
#include "dijkl2.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_dijkl1() {
    std::printf("Test 1: dijkl1 (nati=1, 2 atoms x 1 AO, nmos=2)\n");
    meci_C::nmos = 2;
    molkst_C::norbs = 2;
    molkst_C::numat = 2;
    common_arrays_C::nfirst.assign({0, 1, 2});
    common_arrays_C::nlast.assign({0, 1, 2});
    // c(ip,i) col-major: MO1=(0.6,0.8), MO2=(0.8,-0.6)
    std::vector<std::vector<double>> cv;
    cv.push_back({0.6, 0.8});  // ip=1 -> (i=1,2)
    cv.push_back({0.8, -0.6}); // ip=2
    // w: 2-center nati(atom1)-atom2 integral = 0.3 at w[1]
    std::vector<double> w(8, 0.0);
    w[1] = 0.3;
    std::vector<double> cij(8, 0.0), wcij(8, 0.0), ckl(8, 0.0);
    std::vector<double> xy(16, 0.0);
    dijkl1(cv, 2, 1, w, cij, wcij, ckl, xy);
    // wcij is a workspace reused across (i,j); final state is (i,j)=(2,2):
    // cij atom2=0.72, nati=1.28 -> wcij(1)=wcb=1.28*.3*.25=.096, wcij(2)=wca=.72*.3*.25=.054
    CHK_D(wcij[0], 0.096, 1e-12, "wcij(1) final (atom2 side)");
    CHK_D(wcij[1], 0.054, 1e-12, "wcij(2) final (nati side)");
    // xy(1,1,1,1) = 1.28*.054 + .72*.096 = .13824
    CHK_D(xy[0], 0.13824, 1e-12, "xy(1,1,1,1)");
    // i=2,j=1: xy(2,1,1,1)=.04032, xy(2,1,2,1)=-.13824, xy(2,1,2,2)=-.04032
    CHK_D(xy[1], 0.04032, 1e-12, "xy(2,1,1,1)");
    CHK_D(xy[5], -0.13824, 1e-12, "xy(2,1,2,1)");
    CHK_D(xy[13], -0.04032, 1e-12, "xy(2,1,2,2)");
    // symmetry: xy(1,1,2,2)=xy(2,1,1,2)... use xy(l,k,j,i) mirror: xy(1,2,1,1)? not in (j<=i) loop; check mirror of (2,1,2,1): xy(2,1,1,2)? k,l swapped -> xy(2,1,2,1) already set; check xy(1,2,1,2)=? via (i,j,k,l)=(2,1,2,2) mirror xy(1,2,2,2)
    CHK_D(xy[14], -0.04032, 1e-12, "xy(1,2,2,2) mirror");
}

static void test_dijkl2() {
    std::printf("Test 2: dijkl2 (norbs=2, nmos=2)\n");
    meci_C::nmos = 2;
    molkst_C::norbs = 2;
    // dijkl(k,l,ij) 1-based: size norbs*nmos*3
    meci_C::dijkl.assign(2 * 2 * 3, 0.0);
    double d[12] = {0.13824, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14};
    for (int t = 0; t < 12; ++t) meci_C::dijkl[t] = d[t];
    meci_C::xy.assign(16, 0.0);
    // dc(:,1)={0.1,0.2}, dc(:,2)={0.3,0.4}
    std::vector<std::vector<double>> dc;
    dc.push_back({0.1, 0.2});
    dc.push_back({0.3, 0.4});
    dijkl2(dc);
    CHK_D(meci_C::xy[0], 0.087296, 1e-12, "xy(1,1,1,1) 4x");
    CHK_D(meci_C::xy[1], 0.120472, 1e-12, "xy(2,1,1,1)");
    CHK_D(meci_C::xy[5], 0.164, 1e-12, "xy(2,1,2,1)");
    CHK_D(meci_C::xy[13], 0.256, 1e-12, "xy(2,1,2,2)");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08d batch - dijkl1/dijkl2\n");
    test_dijkl1();
    test_dijkl2();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
