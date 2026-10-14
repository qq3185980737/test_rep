// test_addhb.cpp — unit test for addhb().
//
// addhb only orchestrates hbonds / diagg2 (not yet ported). This test supplies
// stubs that record the calls, so we verify addhb's control flow, the hblims
// threshold selection, and the eigs(nocc1+1:) pointer arithmetic.
//
// Build with:
//   cl /std:c++17 /utf-8 /EHsc test_addhb.cpp addhb.cpp molkst_C.cpp common_arrays_C.cpp
// Exit code 0 == all PASS.

#include <cstdio>
#include <vector>

#include "addhb.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

namespace {

int g_hbonds_called = 0;
double g_cutoff_seen = -1.0;
int g_nij_to_report = 0;
int g_diagg2_called = 0;
const double* g_eigv_seen = nullptr;
int g_memory_error_called = 0;

}  // namespace

// Stub implementations of the not-yet-ported routines.
void memory_error(const char*) { ++g_memory_error_called; }

void hbonds(const double*, int, int, int*, int& nij_loc, double cutoff) {
    ++g_hbonds_called;
    g_cutoff_seen = cutoff;
    nij_loc = g_nij_to_report;
}

void diagg2(int, int, const double* eigv, int*, char*, int, int, double*,
            double*) {
    ++g_diagg2_called;
    g_eigv_seen = eigv;
}

int failures = 0;
void check(const char* name, bool ok) {
    std::printf("%-40s %s\n", name, ok ? "PASS" : "FAIL");
    if (!ok) {
        ++failures;
    }
}

int main() {
    molkst_C::numat = 2;
    molkst_C::norbs = 4;
    // 1-based padded module arrays: [0] is padding.
    common_arrays_C::f.assign(10, 0.0);
    common_arrays_C::eigs = {-999.0, -2.0, -1.0, 3.0, 5.0};  // eigs[1..4]

    int nij = -1;

    // ---- Case 1: nhb=2 -> cutoff 0.1; fake hbonds reports nij=0. ----
    g_nij_to_report = 0;
    addhb(/*nocc1=*/2, /*nvir1=*/2, /*idiagg=*/1, nij, /*nhb=*/2);
    check("case1 hbonds called", g_hbonds_called == 1);
    check("case1 cutoff == hblims(2)=0.1", g_cutoff_seen == 0.1);
    check("case1 diagg2 NOT called (nij==0)", g_diagg2_called == 0);
    check("case1 nij propagated out == 0", nij == 0);
    check("case1 memory_error NOT called", g_memory_error_called == 0);

    // ---- Case 2: nhb=3 -> cutoff 0.01; fake hbonds reports nij=2. ----
    g_nij_to_report = 2;
    addhb(/*nocc1=*/2, /*nvir1=*/2, /*idiagg=*/1, nij, /*nhb=*/3);
    check("case2 hbonds called again", g_hbonds_called == 2);
    check("case2 cutoff == hblims(3)=0.01", g_cutoff_seen == 0.01);
    check("case2 diagg2 called once", g_diagg2_called == 1);
    check("case2 nij propagated out == 2", nij == 2);
    // eigv(1) must alias eigs(nocc1+1) = eigs[3] = 3.0.
    check("case2 eigv(1) aliases eigs[nocc1+1]",
          g_eigv_seen != nullptr && g_eigv_seen[1] == common_arrays_C::eigs[3]);
    check("case2 alias value is 3.0",
          g_eigv_seen != nullptr && g_eigv_seen[1] == 3.0);

    if (failures == 0) {
        std::printf("\nAll checks PASSED.\n");
        return 0;
    }
    std::printf("\n%d check(s) FAILED.\n", failures);
    return 1;
}
