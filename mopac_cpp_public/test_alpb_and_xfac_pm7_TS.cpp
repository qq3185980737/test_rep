// test_alpb_and_xfac_pm7_TS.cpp — verify PM7-TS ALPB table (2056 entries,
// TS-parameterized values differ from PM7 in 58 places).
#include <cmath>
#include <cstdio>
#include "alpb_and_xfac.h"
#include "parameters_C.h"
using namespace parameters_C;
static bool chk(const char* name, double got, double want) {
    bool ok = std::fabs(got - want) < 1e-9;
    std::printf("%-12s got=%g want=%g  %s\n", name, got, want, ok ? "PASS" : "FAIL");
    return ok;
}
int main() {
    alpb_and_xfac_pm7_TS();
    bool ok = true;
    ok &= chk("alpb[1][1]", alpb[1][1], 5.434919);
    ok &= chk("xfac[1][1]", xfac[1][1], 1.449515);
    // a TS-specific difference vs PM7: (xfac,15,7)
    ok &= chk("xfac[15][7]", xfac[15][7], 1.588417);
    ok &= chk("alpb[15][7]", alpb[15][7], 1.594467);
    ok &= chk("alpb[87][87]", alpb[87][87], 1.579660);
    ok &= chk("xfac[87][87]", xfac[87][87], 0.761560);
    ok &= chk("alpb[5][5]", alpb[5][5], 2.181999);
    return ok ? 0 : 1;
}
