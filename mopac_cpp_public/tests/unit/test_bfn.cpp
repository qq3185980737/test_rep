// test_bfn.cpp
#include <cmath>
#include <cstdio>

#include "bfn.h"

int main() {
    double bf[14];
    // x=0: bf(1) = 2*mod(1,2)/1 = 2.0
    bfn(0.0, bf);
    bool r1 = std::fabs(bf[1] - 2.0) < 1e-12;
    std::printf("bf(0) bf[1]=%.4f (expect 2.0) %s\n", bf[1], r1 ? "PASS" : "FAIL");

    // x=4: bf(1) = (e^4 - e^-4)/4
    bfn(4.0, bf);
    double expect = (std::exp(4.0) - std::exp(-4.0)) / 4.0;
    bool r2 = std::fabs(bf[1] - expect) < 1e-9;
    std::printf("bf(4) bf[1]=%.6f (expect %.6f) %s\n", bf[1], expect, r2 ? "PASS" : "FAIL");

    // recurrence consistency: bf(2) = (1*bf(1) - e^4 - e^-4)/4
    double bf2exp = (1.0 * bf[1] - std::exp(4.0) - std::exp(-4.0)) / 4.0;
    bool r3 = std::fabs(bf[2] - bf2exp) < 1e-9;
    std::printf("bf(4) bf[2]=%.6f (expect %.6f) %s\n", bf[2], bf2exp, r3 ? "PASS" : "FAIL");

    return (r1 && r2 && r3) ? 0 : 1;
}
