// test_m03_batchB.cpp — M03 batch B: jab / kab / ijbo numerical comparison
// Expected values computed independently (python, F90 formulas), not from the C++ tables.
#include <cmath>
#include <cstdio>
#include <vector>
#include "jab.h"
#include "kab.h"
#include "ijbo.h"
#include "MOZYME_C.h"
#include "overlaps_C.h"
#include "common_arrays_C.h"

static int passed = 0, failed = 0;
static void check(const char* name, bool ok) {
    if (ok) { ++passed; std::printf("  [PASS] %s\n", name); }
    else { ++failed; std::printf("  [FAIL] %s\n", name); }
}
static bool close(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}

int main() {
    std::printf("M03 batch B tests (jab / kab / ijbo)\n");
    {
        double w[101];
        for (int k = 1; k <= 100; ++k) w[k] = 0.01*k + 0.001*k*k;
        double pja[17], pjb[17];
        for (int c = 1; c <= 16; ++c) { pja[c] = 0.5 + c*0.05; pjb[c] = 1.0 + c*0.03; }
        double f[40] = {0};
        jab(2, 3, pja, pjb, w, f);
        bool ok = true;
        ok = ok && close(f[3], 2.01672);
        ok = ok && close(f[5], 8.39072);
        ok = ok && close(f[6], 78.14752);
        ok = ok && close(f[8], 33.18672);
        ok = ok && close(f[9], 112.70192);
        ok = ok && close(f[10], 136.89592);
        ok = ok && close(f[12], 100.50072);
        ok = ok && close(f[13], 195.60552);
        ok = ok && close(f[14], 231.90672);
        ok = ok && close(f[15], 272.25352);
        ok = ok && close(f[18], 70.1692);
        ok = ok && close(f[19], 72.0732);
        ok = ok && close(f[20], 74.0068);
        ok = ok && close(f[21], 75.97);
        check("jab ia=2,ja=3 accumulation", ok);
    }
    {
        double w[101];
        for (int k = 1; k <= 100; ++k) w[k] = 0.01*k + 0.001*k*k;
        double pk[17];
        for (int c = 1; c <= 16; ++c) pk[c] = 2.0 + c*0.07;
        double f[40] = {0};
        kab(5, 2, pk, w, f);
        bool ok = true;
        ok = ok && close(f[12], -76.16364);
        ok = ok && close(f[13], -79.2166);
        ok = ok && close(f[14], -83.83848);
        ok = ok && close(f[15], -91.98272);
        ok = ok && close(f[17], -110.43004);
        ok = ok && close(f[18], -114.3134);
        ok = ok && close(f[19], -120.16968);
        ok = ok && close(f[20], -130.37312);
        ok = ok && close(f[23], -167.40084);
        ok = ok && close(f[24], -172.485);
        ok = ok && close(f[25], -180.12568);
        ok = ok && close(f[26], -193.30592);
        ok = ok && close(f[30], -303.67964);
        ok = ok && close(f[31], -310.7726);
        ok = ok && close(f[32], -321.39848);
        ok = ok && close(f[33], -339.55872);
        check("kab ia>ja (5,2) subtraction", ok);
    }
    {
        double w[101];
        for (int k = 1; k <= 100; ++k) w[k] = 0.01*k + 0.001*k*k;
        double pk[17];
        for (int c = 1; c <= 16; ++c) pk[c] = 2.0 + c*0.07;
        double f[25] = {0};
        kab(2, 3, pk, w, f);
        bool ok = true;
        ok = ok && close(f[5], -76.16364);
        ok = ok && close(f[6], -110.43004);
        ok = ok && close(f[7], -167.40084);
        ok = ok && close(f[8], -382.89624);
        ok = ok && close(f[9], -114.3134);
        ok = ok && close(f[10], -172.485);
        ok = ok && close(f[11], -310.7726);
        ok = ok && close(f[12], -83.83848);
        ok = ok && close(f[13], -120.16968);
        ok = ok && close(f[14], -180.12568);
        ok = ok && close(f[15], -321.39848);
        ok = ok && close(f[17], -91.98272);
        ok = ok && close(f[18], -130.37312);
        ok = ok && close(f[19], -193.30592);
        ok = ok && close(f[20], -339.55872);
        check("kab ia<=ja (2,3) triangle", ok);
    }
    {
        // lijbo fast path
        MOZYME_C::lijbo = true;
        MOZYME_C::nijbo.assign(4, std::vector<int>(4, 0));
        MOZYME_C::nijbo[2][3] = 5;
        bool ok = (ijbo(2, 3) == 5);
        MOZYME_C::lijbo = false;
        // lookup table for atom pair table: iij/numij 1-based, ijall/iijj 1-based (0 padding)
        MOZYME_C::iij.assign(4, 0); MOZYME_C::numij.assign(4, 0);
        MOZYME_C::ijall.assign(4, 0); MOZYME_C::iijj.assign(4, 0);
        MOZYME_C::iij[2] = 1; MOZYME_C::numij[2] = 2;   // atom 2: 2 partners
        MOZYME_C::ijall[1] = 1; MOZYME_C::ijall[2] = 2; // partners (1,2)
        MOZYME_C::iijj[1] = 10; MOZYME_C::iijj[2] = 20;
        // close pair (r <= cutof2) -> binary search hits
        common_arrays_C::coord.assign(4, std::vector<double>(4, 0.0));
        overlaps_C::cutof1 = 4.0; overlaps_C::cutof2 = 2.0;
        ok = ok && (ijbo(2, 1) == 10);  // ind_i=2, ind_j=1 -> iijj[1]
        ok = ok && (ijbo(1, 2) == 10);  // ind_j is always min: ind_i=2, ind_j=1 -> iijj[1]
        ok = ok && (ijbo(2, 2) == 20);  // ind_i=2, ind_j=2 -> iijj[2]
        // medium distance: cutof2 < r <= cutof1 -> -2
        common_arrays_C::coord[1][2] = 1.5;
        ok = ok && (ijbo(1, 2) == -2);
        // far: r > cutof1 -> -1
        common_arrays_C::coord[1][2] = 3.0;
        ok = ok && (ijbo(1, 2) == -1);
        check("ijbo fast/lookup/distance paths", ok);
    }

    std::printf("M03 batch B: %d passed, %d failed\n", passed, failed);
    return failed ? 1 : 0;
}
