// test_resolv.cpp — verify degenerate-LMO resolution.
// 1) Single-s-orbital atoms (nlast==nfirst) never trigger the candidate path:
//    c must come back unchanged.
// 2/3) Similar LMOs on multi-orbital atoms form degenerate pairs: resolv
//    rotates them so the off-diagonal F interaction <i|F|j> vanishes
//    (diagonalizes the small secular matrix) while columns stay in the
//    original span (transform is invertible).
#include <cstdio>
#include <cmath>
#include <vector>
#include "resolv.h"
#include "rsp.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using molkst_C::norbs;
using molkst_C::numat;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    if (std::fabs(got - want) > tol) {
        ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.10g want %.10g\n", name, got, want);
    }
}

int main() {
    // Case 1: H2-like, single s orbital per atom -> no candidate -> unchanged.
    {
        numat = 2; norbs = 2;
        nfirst = {0, 1, 2}; nlast = {0, 1, 2};
        double c[4] = {0.9, 0.4, 0.2, 0.8};
        double cold[4] = {0.9, 0.4, 0.2, 0.8};
        double eig[2] = {-0.5, -0.3};
        double c0[4];
        for (int i = 0; i < 4; ++i) c0[i] = c[i];
        resolv(c, cold, 2, eig, 2);
        bool same = true;
        for (int i = 0; i < 4; ++i) if (c[i] != c0[i]) same = false;
        chk(same, "single-orbital atoms leave c unchanged");
    }
    // Case 2: two similar LMOs (0.8,0.6)/(0.6,0.8) -> F cross term -> 0.
    {
        numat = 1; norbs = 2;
        nfirst = {0, 1}; nlast = {0, 2};
        double cold[4] = {1, 0, 0, 1};
        double eig[2] = {1.0, 2.0};
        double c[4] = {0.8, 0.6, 0.6, 0.8};
        resolv(c, cold, 2, eig, 2);
        double cross = c[0]*1.0*c[2] + c[1]*2.0*c[3];
        chk_rel(cross, 0.0, 1e-10, "2x2 F cross term zero");
        double det = c[0]*c[3] - c[1]*c[2];
        chk(std::fabs(det) > 1e-6, "2x2 rotation invertible (span preserved)");
        double d1 = c[0]*1.0*c[0] + c[1]*2.0*c[1];
        double d2 = c[2]*1.0*c[2] + c[3]*2.0*c[3];
        chk_rel(d1, 0.053210, 1e-4, "2x2 diag1 = small eigenvalue");
        chk_rel(d2, 2.946790, 1e-4, "2x2 diag2 = large eigenvalue");
    }
    // Case 3: three similar LMOs -> 3x3, all F cross terms zero, diagonal
    // terms pairwise distinct (diagonalized).
    {
        numat = 1; norbs = 3;
        nfirst = {0, 1}; nlast = {0, 3};
        double cold[9] = {1,0,0, 0,1,0, 0,0,1};
        double eig[3] = {1.0, 2.0, 3.0};
        double c[9] = {0.8,0.6,0.0, 0.0,0.8,0.6, 0.6,0.0,0.8};
        resolv(c, cold, 3, eig, 3);
        double cross12 = c[0]*1.0*c[3] + c[1]*2.0*c[4] + c[2]*3.0*c[5];
        double cross13 = c[0]*1.0*c[6] + c[1]*2.0*c[7] + c[2]*3.0*c[8];
        double cross23 = c[3]*1.0*c[6] + c[4]*2.0*c[7] + c[5]*3.0*c[8];
        chk_rel(cross12, 0.0, 1e-10, "3x3 F cross 1-2 zero");
        chk_rel(cross13, 0.0, 1e-10, "3x3 F cross 1-3 zero");
        chk_rel(cross23, 0.0, 1e-10, "3x3 F cross 2-3 zero");
        double d[3];
        for (int j = 0; j < 3; ++j)
            d[j] = c[j*3]*1.0*c[j*3] + c[j*3+1]*2.0*c[j*3+1] + c[j*3+2]*3.0*c[j*3+2];
        bool distinct = true;
        for (int a = 0; a < 3; ++a)
            for (int b = a + 1; b < 3; ++b)
                if (std::fabs(d[a] - d[b]) < 1e-6) distinct = false;
        chk(distinct, "3x3 diagonal terms distinct (diagonalized)");
    }
    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}