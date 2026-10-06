// test_switch.cpp — verify switch_method() loads each Hamiltonian's reference
// parameter set into parameters_C and applies the ALPB / fractional-metal /
// SPARKLES post-processing.
#include "switch.h"

#include "molkst_C.h"
#include "parameters_C.h"
#include "parameters_for_PM7_Sparkles_C.h"

#include <cstdio>

using namespace parameters_C;
using namespace molkst_C;
using namespace parameters_for_PM7_Sparkles_C;

static bool near(double a, double b) {
    double d = a - b; if (d < 0) d = -d;
    double s = a < 0 ? -a : a; if (b < 0) s = -s;
    return d <= 1e-8 || d <= 1e-8 * s;
}

static void reset_flags() {
    method_mndo = method_pm3 = method_pm6 = method_pm7 = method_pm7_ts =
        method_mndod = method_rm1 = false;
    keywrd.clear();
}

int main() {
    bool ok = true;
    auto chk = [&](bool cond, const char* msg) {
        if (!cond) { std::printf("FAIL %s\n", msg); ok = false; }
    };

    // ---- AM1 (default branch) ----
    reset_flags();
    switch_method();
    chk(near(uss[1], -11.3964270), "AM1 uss[1]");
    chk(near(betas[1], -6.1737870), "AM1 betas[1]");
    chk(near(zs[1], 1.1880780), "AM1 zs[1]");
    chk(near(gss[1], 12.8480000), "AM1 gss[1]");
    chk(near(alp[1], 2.8823240), "AM1 alp[1]");
    chk(near(polvol[1], 0.1719890), "AM1 polvol[1]");
    chk(near(guess1[1][1], 0.1227960), "AM1 guess1[1][1]");
    chk(near(pocord[42], 1.334), "AM1 pocord[42] Mo");
    chk(near(alpb[3][1], 2.9751160), "AM1 alpb[3][1]");
    chk(near(alpb[1][3], 2.9751160), "AM1 alpb symmetrized");
    chk(near(tore[85], -0.5) && near(tore[87], 0.5), "fractional metal tore");
    chk(near(alp[85], 3.0) && near(gss[85], 10.0), "fractional metal params");
    chk(near(f0sd_store[1], f0sd[1]) && near(g2sd_store[1], g2sd[1]), "f0sd/g2sd store");

    // ---- PM3 ----
    reset_flags(); method_pm3 = true;
    switch_method();
    chk(near(uss[1], -13.0733210), "PM3 uss[1]");

    // ---- MNDO ----
    reset_flags(); method_mndo = true;
    switch_method();
    chk(near(uss[1], -11.9062760), "MNDO uss[1]");

    // ---- MNDO/d ----
    reset_flags(); method_mndod = true;
    switch_method();
    chk(near(uss[1], -11.9062760), "MNDOD uss[1]");
    chk(near(guess1[6][1], 0.0), "MNDOD guess1 zeroed");
    chk(near(f0sd[1], 0.0), "MNDOD f0sd zeroed");

    // ---- PM6 ----
    reset_flags(); method_pm6 = true;
    switch_method();
    chk(near(uss[1], -11.2469580), "PM6 uss[1]");
    chk(near(alpb[1][1], 3.5409420), "PM6 alpb[1][1]");
    chk(near(xfac[1][1], 2.2435870), "PM6 xfac[1][1]");

    // ---- PM7 ----
    reset_flags(); method_pm7 = true;
    switch_method();
    chk(near(uss[1], -11.0701120), "PM7 uss[1]");
    chk(near(alp[13], 5.3416850), "PM7 alp[13] (first alp7 data)");
    chk(near(alpb[1][1], 4.0511630), "PM7 alpb[1][1]");
    chk(near(xfac[1][1], 2.8456270), "PM7 xfac[1][1]");

    // ---- PM7 SPARKLES ----
    reset_flags(); method_pm7 = true; keywrd = "PM7 SPARK";
    switch_method();
    chk(near(alp[57], alp7sp[57]), "PM7 sparkle alp[57]");
    chk(near(gss[57], gss7sp[57]), "PM7 sparkle gss[57]");
    chk(near(uss[57], 0.0) && near(zs[57], 0.0) && near(zp[57], 0.0), "PM7 sparkle zeroed uss/zs/zp");
    chk(near(guess1[57][1], gues7sp1[57][1]), "PM7 sparkle guess1[57][1]");
    chk(near(guess1[57][3], 0.0), "PM7 sparkle guess1[57][3] zeroed");
    chk(near(alpb[57][1], 0.0) && near(xfac[1][71], 0.0), "PM7 sparkle alpb/xfac lanthanide zeroed");

    // ---- PM7-TS ----
    reset_flags(); method_pm7_ts = true;
    switch_method();
    chk(near(uss[1], -11.2617750), "PM7_TS uss[1]");

    // ---- RM1 ----
    reset_flags(); method_rm1 = true;
    switch_method();
    chk(near(uss[1], -11.9606770), "RM1 uss[1]");
    chk(near(zs[1], 1.0826737), "RM1 zs[1]");

    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
