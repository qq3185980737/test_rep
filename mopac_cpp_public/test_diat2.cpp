// test_diat2.cpp — numerical verification of the full diat.F90 translation
// (diat driver + diat2 sp-path + ss general STO overlap). Reference values
// computed independently in gen_test_diat.py from the F90 formulas.
#include <cmath>
#include <cstdio>
#include <vector>

#include "diat.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"

using namespace funcon_C;
using namespace molkst_C;
using namespace overlaps_C;
using namespace parameters_C;

static int failures = 0;
static void check(const char* name, double got, double want, double tol = 1e-8) {
    bool ok = std::fabs(got - want) <= tol * std::max(1.0, std::fabs(want));
    std::printf("%-16s = %16.10f  want %16.10f  %s\n", name, got, want,
                ok ? "PASS" : "FAIL");
    if (!ok) ++failures;
}

int main() {
    numcal = 1;
    a0 = 1.0;          // simplify: a0 = 1 bohr
    cutof1 = 1e6;      // no early cutoff
    cutof2 = 1e-4;
    // Element data: H (1) and C (6).
    natorb[1] = 1; zs[1] = 1.0; zp[1] = 0.0; zd[1] = 0.0;
    natorb[6] = 4; zs[6] = 2.0; zp[6] = 1.5; zd[6] = 0.0;
    npq[1][1] = 1; npq[1][2] = 1; npq[1][3] = 0;
    npq[6][1] = 2; npq[6][2] = 2; npq[6][3] = 0;

    // --- ss (general STO overlap) ---
    check("ss(1s,1s)", ss_overlap(1, 1, 1, 1, 1, 1.0, 1.0, 1.0, 1.0), 0.858385362733);
    check("ss(2s,2s)", ss_overlap(2, 2, 1, 1, 1, 1.0, 1.0, 1.0, 1.0), 0.948311448353);
    check("ss(2ps,2ps)", ss_overlap(2, 2, 2, 2, 1, 1.0, 1.0, 1.0, 1.0), -0.735758882343);
    check("ss(s,s,m=1)=0", ss_overlap(1, 1, 1, 1, 2, 1.0, 1.0, 1.0, 1.0), 0.0);
    check("ss(2pp,2pp)", ss_overlap(2, 2, 2, 2, 2, 1.0, 1.0, 1.0, 1.0), 0.907435954890);

    // --- diat: H-H (diat2 case 1), along x at 1 A ---
    std::vector<std::vector<double>> di(10, std::vector<double>(10, 0.0));
    diat(1, 1, std::vector<double>{1.0, 0.0, 0.0}, di);
    check("diHH(1,1)", di[1][1], 0.858385362733);
    double ssum = 0.0;
    for (int i = 1; i <= 9; ++i)
        for (int j = 1; j <= 9; ++j) ssum += std::fabs(di[i][j]);
    check("diHH total", ssum, 0.858385362733);
    // di is symmetric: di(1,1) only nonzero for H-H along x
    check("diHH(2,2)=0", di[2][2], 0.0);

    // --- diat: H-C (diat2 case 2), along x at 1 A ---
    for (auto& row : di) std::fill(row.begin(), row.end(), 0.0);
    diat(1, 6, std::vector<double>{1.0, 0.0, 0.0}, di);
    check("diHC(1,1)", di[1][1], 0.815351539118);
    check("diHC(1,3)=0", di[1][3], 0.0);  // pi component projects to 0 along x

    // --- diat: C-C (diat2 case 4), along x at 1 A ---
    for (auto& row : di) std::fill(row.begin(), row.end(), 0.0);
    diat(6, 6, std::vector<double>{1.0, 0.0, 0.0}, di);
    check("diCC(1,1)", di[1][1], 0.815019150158);
    ssum = 0.0;
    for (int i = 1; i <= 9; ++i)
        for (int j = 1; j <= 9; ++j) ssum += std::fabs(di[i][j]);
    check("diCC total", ssum, 3.930275167320);

    // --- cutoff behavior: far apart -> all zeros ---
    for (auto& row : di) std::fill(row.begin(), row.end(), 0.0);
    diat(1, 1, std::vector<double>{100.0, 0.0, 0.0}, di);
    ssum = 0.0;
    for (int i = 1; i <= 9; ++i)
        for (int j = 1; j <= 9; ++j) ssum += std::fabs(di[i][j]);
    check("diHH r=100 -> 0", ssum, 0.0);

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASS",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
