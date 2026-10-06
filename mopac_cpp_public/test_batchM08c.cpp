// test_batchM08c.cpp  M08 ijkl end-to-end (cij/ckl + partxy + dijkl + xy)
#include <cstdio>
#include <cmath>
#include <vector>
namespace cosmo_C { extern bool iseps; }
#include "ijkl.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08c batch - ijkl\n");
    // 2 atoms x 1 AO each; nmos=1, nelec=0
    molkst_C::norbs = 2;
    molkst_C::numat = 2;
    molkst_C::lm61 = 8;
    common_arrays_C::nfirst.assign({0, 1, 2});
    common_arrays_C::nlast.assign({0, 1, 2});
    // w: w[1]=0.4 atom1 1-center, w[2]=0.3 2-center, w[3]=0.2 atom2 1-center
    common_arrays_C::w.assign(16, 0.0);
    common_arrays_C::w[1] = 0.4; common_arrays_C::w[2] = 0.3; common_arrays_C::w[3] = 0.2;
    extern bool cosmo_C_iseps();
    cosmo_C::iseps = false;

    // cp = active MO (norbs x nmos col-major): cp(1,1)=0.6, cp(2,1)=0.8
    double cp[2] = {0.6, 0.8};
    // cf = all MOs (norbs x norbs): col0 = cp, col1 = 0.8,-0.6
    double cf[4] = {0.6, 0.8, 0.8, -0.6};
    int mpair = 1;              // nmos*(nmos+1)/2 = 1
    std::vector<double> dijkl(2 * 1 * mpair, 0.0);
    double cij[8] = {}, ckl[8] = {}, wcij[8] = {};
    double xy[1] = {0.0};
    ijkl(cp, cf, 0, 1, dijkl.data(), cij, ckl, wcij, xy);

    // Hand-checked: wcij(1)=cij1*w1*.25 + cij2*w2*.25 = .072+.096=.168
    //             wcij(2)=cij2*w3*.25 + cij1*w2*.25 = .064+.054=.118
    CHK_D(wcij[0], 0.168, 1e-12, "wcij(1) via partxy");
    CHK_D(wcij[1], 0.118, 1e-12, "wcij(2) via partxy");
    // dijkl(1,1,1) = .72*.168 + 1.28*.118 = .272
    CHK_D(dijkl[0], 0.272, 1e-12, "dijkl(1,1,1)");
    // dijkl(2,1,1) = .96*.168 - .96*.118 = .048
    CHK_D(dijkl[1], 0.048, 1e-12, "dijkl(2,1,1)");
    // xy(1,1,1,1) = dijkl(1,1,1)
    CHK_D(xy[0], 0.272, 1e-12, "xy(1,1,1,1) spread");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
