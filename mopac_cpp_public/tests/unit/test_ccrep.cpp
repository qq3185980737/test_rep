// test_ccrep.cpp
#include <cmath>
#include <cstdio>

#include "ccrep.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

int main() {
    using namespace parameters_C;
    funcon_C::a0 = 0.529177;
    molkst_C::method_pm6 = true;
    tore[6] = 4.0;
    xfac[6][6] = 0.0;
    alpb[6][6] = 0.0;

    double r = 10.0;  // Bohr, large enough to suppress LJ
    double enuclr = 0.0;
    ccrep(6, 6, r, 1.0, enuclr);

    double expected = std::abs(10.0 * exp(-2.18 * (10.0 * funcon_C::a0)) * 16.0) + 16.0;
    bool ok = std::fabs(enuclr - expected) < 1e-6;
    std::printf("ccrep=%.4f expected~%.4f %s\n", enuclr, expected,
                ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
