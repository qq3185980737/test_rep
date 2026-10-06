// test_dcart.cpp
#include <cmath>
#include <cstdio>
#include "dcart.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal = 1; clower = 3.0; cupper = 5.0; cutofp = 4.0;
    double d1 = derp(1.0);     // < clower -> -1
    double d2 = derp(10.0);    // > cupper -> 0
    bool ok = (std::abs(d1 + 1.0) < 1e-9 && d2 == 0.0);
    std::printf("derp(1)=%g derp(10)=%g %s\n", d1, d2, ok?"PASS":"FAIL");
    return ok?0:1;
}
