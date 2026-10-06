// test_params_batch.cpp - spot-check first data value of each parameters_for_* module.
#include <cstdio>
#include "parameters_for_AM1_C.h"
#include "parameters_for_AM1_Sparkles_C.h"
#include "parameters_for_mndo_C.h"
#include "parameters_for_mndod_C.h"
#include "parameters_for_PM3_C.h"
#include "parameters_for_PM3_Sparkles_C.h"
#include "parameters_for_PM6_C.h"
#include "parameters_for_PM6_Sparkles_C.h"
#include "parameters_for_PM7_C.h"
#include "parameters_for_PM7_Sparkles_C.h"
#include "parameters_for_PM7_TS_C.h"
#include "parameters_for_RM1_C.h"
#include "parameters_for_RM1_Sparkles_C.h"

using namespace parameters_for_AM1_C;
using namespace parameters_for_AM1_Sparkles_C;
using namespace parameters_for_mndo_C;
using namespace parameters_for_mndod_C;
using namespace parameters_for_PM3_C;
using namespace parameters_for_PM3_Sparkles_C;
using namespace parameters_for_PM6_C;
using namespace parameters_for_PM6_Sparkles_C;
using namespace parameters_for_PM7_C;
using namespace parameters_for_PM7_Sparkles_C;
using namespace parameters_for_PM7_TS_C;
using namespace parameters_for_RM1_C;
using namespace parameters_for_RM1_Sparkles_C;

int main() {
    bool ok = true;
    auto near = [](double a, double b) { return (a - b) * (a - b) < 1e-8; };
    if (!near(ussam1[1], -11.3964270e0)) { std::printf("FAIL AM1: ussam1[1] = %g\n", ussam1[1]); ok = false; }
    if (!near(alpam1sp[57], 2.1879021e0)) { std::printf("FAIL AM1_Sparkles: alpam1sp[57] = %g\n", alpam1sp[57]); ok = false; }
    if (!near(ussm[1], -11.9062760e0)) { std::printf("FAIL mndo: ussm[1] = %g\n", ussm[1]); ok = false; }
    if (!near(ussd[1], -11.9062760e0)) { std::printf("FAIL mndod: ussd[1] = %g\n", ussd[1]); ok = false; }
    if (!near(usspm3[1], -13.0733210e0)) { std::printf("FAIL PM3: usspm3[1] = %g\n", usspm3[1]); ok = false; }
    if (!near(alpPM3sp[57], 2.0790327441e0)) { std::printf("FAIL PM3_Sparkles: alpPM3sp[57] = %g\n", alpPM3sp[57]); ok = false; }
    if (!near(uss6[1], -11.246958e0)) { std::printf("FAIL PM6: uss6[1] = %g\n", uss6[1]); ok = false; }
    if (!near(alp6sp[57], 2.0955474333e0)) { std::printf("FAIL PM6_Sparkles: alp6sp[57] = %g\n", alp6sp[57]); ok = false; }
    if (!near(uss7[1], -11.070112e0)) { std::printf("FAIL PM7: uss7[1] = %g\n", uss7[1]); ok = false; }
    if (!near(alp7sp[57], 2.72878366e0)) { std::printf("FAIL PM7_Sparkles: alp7sp[57] = %g\n", alp7sp[57]); ok = false; }
    if (!near(uss7_TS[1], -11.261775e0)) { std::printf("FAIL PM7_TS: uss7_TS[1] = %g\n", uss7_TS[1]); ok = false; }
    if (!near(ussRM1[1], -11.9606770e0)) { std::printf("FAIL RM1: ussRM1[1] = %g\n", ussRM1[1]); ok = false; }
    if (!near(alprm1sp[57], 2.06878895e0)) { std::printf("FAIL RM1_Sparkles: alprm1sp[57] = %g\n", alprm1sp[57]); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
