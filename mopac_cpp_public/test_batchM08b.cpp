// test_batchM08b.cpp  M08 partxy (MNDO two-electron transform) numerics
#include <cstdio>
#include <cmath>
#include <vector>
#include "partxy.h"
#include "molkst_C.h"
#include "common_arrays_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg) do { if (c) { ++n_pass; std::printf("  [PASS] %s\n", msg); } \
  else { ++n_fail; std::printf("  [FAIL] %s\n", msg); } } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a)-(b)) < (tol), msg)

static void test_partxy_2atoms() {
    std::printf("Test 1: partxy 2 atoms x 1 AO each\n");
    molkst_C::numcal = 1;
    molkst_C::numat = 2;
    molkst_C::lm61 = 8;
    common_arrays_C::nfirst.assign({0, 1, 2});
    common_arrays_C::nlast.assign({0, 1, 2});
    // c34: atom1 dist=0.5, atom2 dist=1.5 (1-based -> c34[1], c34[2])
    double c34[8] = {0.0, 0.5, 1.5, 0, 0, 0, 0, 0};
    double pq34[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    // w: w[1]=0.4 (atom1 1-center), w[2]=0.3 (2-center block), w[3]=0.2 (atom2 1-center)
    double w[16] = {0.0, 0.4, 0.3, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    partxy(c34, pq34, w);
    // atom1 1-center: pq34[1] += c34[1]*w[1]*0.25 = 0.5*0.4*0.25 = 0.05
    // 2-center formxy: pq34[2] += cb(0.5)*w[2]*0.25 = 0.0375
    //                 pq34[1] += ca(1.5)*w[2]*0.25 = 0.1125
    // atom2 1-center: pq34[2] += c34[2]*w[3]*0.25 = 1.5*0.2*0.25 = 0.075
    CHK_D(pq34[0], 0.05 + 0.1125, 1e-12, "pq34(1) total");
    CHK_D(pq34[1], 0.0375 + 0.075, 1e-12, "pq34(2) total");
}

static void test_partxy_single_atom() {
    std::printf("Test 2: partxy single atom (1s)\n");
    molkst_C::numcal = 2;
    molkst_C::numat = 1;
    molkst_C::lm61 = 8;
    common_arrays_C::nfirst.assign({0, 1});
    common_arrays_C::nlast.assign({0, 1});
    double c34[8] = {0.0, 2.0, 0, 0, 0, 0, 0, 0};
    double pq34[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double w[16] = {0.0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    partxy(c34, pq34, w);
    CHK_D(pq34[0], 2.0 * 0.5 * 0.25, 1e-12, "single atom pq34(1)");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M08b batch - partxy\n");
    test_partxy_2atoms();
    test_partxy_single_atom();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
