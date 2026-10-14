// test_dang.cpp
#include <cmath>
#include <cstdio>
#include "dang.h"
int main() {
    double a1=1,a2=0,b1=0,b2=1,rc;
    dang(a1,a2,b1,b2,rc);
    bool ok = (std::abs(rc + 4.71238898) < 1e-6);
    std::printf("rcos=%.4f (expect -4.7124) %s\n", rc, ok?"PASS":"FAIL");
    return ok?0:1;
}
