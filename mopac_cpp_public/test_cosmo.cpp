// test_cosmo.cpp
#include <cmath>
#include <cstdio>
#include "cosmo.h"
int main() {
    double aar,abr,ara,arad,arb,arbd,rinc;
    ansude(1.0,1.0,1.0,1.0, aar,abr,ara,arad,arb,arbd,rinc);
    bool ok = (std::abs(ara-9.371135)<1e-3 && std::abs(aar-1.517)<1e-2);
    std::printf("ara=%.6f aar=%.3f %s\n", ara, aar, ok?"PASS":"FAIL");
    return ok?0:1;
}
