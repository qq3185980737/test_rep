// test_alpb_and_xfac_pm6.cpp — verify PM6 ALPB table (1866 entries).
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
    alpb_and_xfac_pm6();
    bool ok = true;
    ok &= chk("alpb[1][1]", alpb[1][1], 3.540942);
    ok &= chk("xfac[1][1]", xfac[1][1], 2.243587);
    ok &= chk("alpb[2][1]", alpb[2][1], 2.989881);
    ok &= chk("xfac[2][1]", xfac[2][1], 2.371199);
    ok &= chk("alpb[2][2]", alpb[2][2], 3.783559);
    ok &= chk("alpb[87][87]", alpb[87][87], 1.579660);
    ok &= chk("xfac[87][87]", xfac[87][87], 0.761560);
    ok &= chk("alpb[5][5]", alpb[5][5], 3.318624);
    return ok ? 0 : 1;
}
