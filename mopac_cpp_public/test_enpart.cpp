#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
// test_enpart.cpp — verification for enpart.cpp (MOPAC 2016 energy
// partitioning).  Scenario: H2, one s orbital per atom.
//   - resonance: e(2,1) = 2*p(2)*h(2) = -23.0 (checked via RESONANCE line;
//     packed lower triangle: p[2]=P(1,2), h[2]=H(1,2))
//   - e(n,2)=g from rotate(); with the current placeholder kernels enuc=0,
//     so N-N is zero; e-e repulsion and exchange come from the same w2
//     stream (0), so the totals are dominated by resonance + one-center.
//   - run completes and prints the SUMMARY block.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "enpart.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace parameters_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}

int main() {
    numat = 2;
    nfirst = {0, 1, 2};
    nlast = {0, 1, 2};
    nat = {0, 1, 1};
    uhf = false;
    numcal = 1;
    itemp_1 = 0;
    keywrd = "";
    p = {0, 1.0, 1.0, 0.25};          // p[1]=P(1,1), p[2]=P(1,2), p[3]=P(2,2)
    h = {0, -11.5, -11.5, -2.0};      // h[1]=H(1,1), h[2]=H(1,2), h[3]=H(2,2)
    pa = p;
    pb = p;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][2] = 0.0; coord[2][2] = 0.0; coord[3][2] = 0.74;
    // post_scf_corrections is called at the end of enpart; use the pm7_minus
    // arm so correction=0 and E_hb=E_disp=0.
    method_pm6_d3h4x = false; method_pm6_d3h4 = false; method_pm6_d3 = false;
    method_pm6_dh_plus = false; method_pm6_dh2 = false; method_pm6_dh2x = false;
    method_pm7_hh = false; method_pm7 = false; method_pm7_minus = true;

    std::freopen("_enpart_out.txt", "w", stdout);
    enpart();
    std::fflush(stdout);
    std::freopen("CON", "w", stdout);

    std::string out;
    {   FILE* f = std::fopen("_enpart_out.txt", "rb");
        std::fseek(f, 0, SEEK_END); long sz = std::ftell(f); std::fseek(f, 0, SEEK_SET);
        out.assign((size_t)sz, '\0'); std::fread(&out[0], 1, (size_t)sz, f); std::fclose(f); }

    // Resonance energy line: " RESONANCE ENERGY" followed by a number = -1.0
    std::size_t pos = out.find("RESONANCE ENERGY");
    chk(pos != std::string::npos, "enpart prints RESONANCE ENERGY");
    if (pos != std::string::npos) {
        std::string tail = out.substr(pos + 16, 20);
        double val = 0.0;
        if (std::sscanf(tail.c_str(), "%lf", &val) == 1)
            chk(std::fabs(val - (-23.0)) < 1e-9, "enpart resonance = 2*p(2)*h(2) = -23");
        else chk(false, "enpart resonance parse");
    }
    chk(out.find("***  SUMMARY OF ENERGY PARTITION  ***") != std::string::npos, "enpart prints summary");
    chk(out.find("ETOT (EONE + ETWO)") != std::string::npos, "enpart prints ETOT");
    chk(out.find("For more detail") != std::string::npos, "enpart non-LARGE hint");
    chk(E_hb == 0.0 && E_disp == 0.0, "enpart post_scf arm resets E_hb/E_disp");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
