// test_alpb_and_xfac_mndo.cpp
#include <cmath>
#include <cstdio>
#include "alpb_and_xfac_mndo.h"
#include "parameters_C.h"

using namespace parameters_C;

static bool chk(const char* name, double got, double want) {
    bool ok = std::fabs(got - want) < 1e-9;
    std::printf("%-12s got=%g want=%g  %s\n", name, got, want, ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    alpb_and_xfac_mndo();
    bool ok = true;
    ok &= chk("alpb[3][1]", alpb[3][1], 1.574938);
    ok &= chk("xfac[3][1]", xfac[3][1], 1.267656);
    ok &= chk("alpb[78][78]", alpb[78][78], 0.999988);
    ok &= chk("xfac[78][78]", xfac[78][78], 0.099990);
    ok &= chk("alpb[53][23]", alpb[53][23], 3.779472);
    ok &= chk("xfac[29][29]", xfac[29][29], 0.129793);
    ok &= chk("alpb[10][10]", alpb[10][10], 0.0);  // unassigned -> 0
    return ok ? 0 : 1;
}
