// test_alpb_and_xfac_am1.cpp
#include <cmath>
#include <cstdio>
#include "alpb_and_xfac_am1.h"
#include "parameters_C.h"

using namespace parameters_C;

static bool chk(const char* name, double got, double want) {
    bool ok = std::fabs(got - want) < 1e-9;
    std::printf("%-12s got=%g want=%g  %s\n", name, got, want, ok ? "PASS" : "FAIL");
    return ok;
}

int main() {
    alpb_and_xfac_am1();
    bool ok = true;
    ok &= chk("alpb[3][1]", alpb[3][1], 2.975116);
    ok &= chk("xfac[3][1]", xfac[3][1], 10.000006);
    ok &= chk("alpb[56][56]", alpb[56][56], 1.110416);
    ok &= chk("xfac[56][56]", xfac[56][56], 0.101890);
    ok &= chk("alpb[42][42]", alpb[42][42], 2.550000);
    ok &= chk("xfac[42][42]", xfac[42][42], 6.000000);
    // an unassigned entry must be zero
    ok &= chk("alpb[2][2]", alpb[2][2], 0.0);
    ok &= chk("xfac[10][10]", xfac[10][10], 0.0);
    return ok ? 0 : 1;
}
