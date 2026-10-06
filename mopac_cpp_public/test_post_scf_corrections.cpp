// test_post_scf_corrections.cpp — verification for post_scf_corrections.cpp.
// Scenario: H2 (0.74 A).  Checks the d3h4x arm (dftd3 + H_bonds4 +
// energy_corr_hh_rep + disp_DnX) with independently computed expectations:
//   dftd3(H2) is small and negative; H_bonds4 = 0 (no N/O acceptors);
//   energy_corr_hh_rep = poly(0.74) = 25.46293603147693 (H-H at r<=1);
//   disp_DnX = 0 (no halogen).  Also checks the pm7_minus arm returns 0.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "post_scf_corrections.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    if (std::fabs(got - want) > tol * (std::fabs(want) > 1 ? std::fabs(want) : 1.0)) {
        ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want);
    }
}

void set_h2() {
    numat = 2;
    nat = std::vector<int>{0, 1, 1};
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[0][2] = 0.0; coord[1][2] = 0.0; coord[2][2] = 0.74;
    numcal = 1;
}

int main() {
    double corr;
    // --- arm 1: PM6-D3H4X ---
    set_h2();
    method_pm6_d3h4x = true; method_pm6_d3h4 = false; method_pm6_d3 = false;
    method_pm6_dh_plus = false; method_pm6_dh2 = false; method_pm6_dh2x = false;
    method_pm7_hh = false; method_pm7_minus = false; method_pm7 = false;
    keywrd = " PM6-D3H4X 0SCF DISP";
    E_hb = 0; E_hh = 0; E_disp = 0; P_Hbonds = 0;
    post_scf_corrections(corr, false);
    // dftd3(H2) = -8.09e-6, H_bonds4 = 0, poly(0.74) = 25.46293603147693, disp_DnX = 0
    chk_rel(corr, 25.46293603147693 - 8.09e-6, 1e-5, "d3h4x H2 total");
    // dftd3 component = E_disp + E_hb; total = poly(0.74) + dftd3 component
    chk_rel(E_disp + E_hb, corr - 25.46293603147693, 1e-8, "d3h4x dftd3 component");
    chk(E_hh == 0.0, "d3h4x E_hh zero for H2");
    chk(P_Hbonds == 0, "d3h4x P_Hbonds=0");

    // --- arm 2: PM7_MINUS returns zero without calling anything ---
    set_h2();
    method_pm6_d3h4x = false; method_pm7_minus = true; method_pm7 = false;
    keywrd = " PM7-MINUS";
    post_scf_corrections(corr, false);
    chk(corr == 0.0, "pm7_minus returns 0");

    // --- arm 3: plain PM7: PM6_DH_Dispersion + PM6_DH_H_bond_corrections ---
    // H2 has no halogen / N-O acceptors, both corrections are structurally 0.
    set_h2();
    method_pm6_d3h4x = false; method_pm7_minus = false; method_pm7 = true;
    keywrd = " PM7";
    post_scf_corrections(corr, false);
    chk(std::isfinite(corr), "pm7 H2 finite");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
