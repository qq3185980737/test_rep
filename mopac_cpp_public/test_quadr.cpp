// test_quadr.cpp
#include <cstdio>
#include <cmath>
#include "quadr.h"
int main() {
    double a,b,c; quadr(0,1,4,1,2,a,b,c);
    bool ok=(std::fabs(a)<1e-9 && std::fabs(b)<1e-9 && std::fabs(c-1)<1e-9);
    std::printf("a=%g b=%g c=%g %s\n",a,b,c,ok?"PASS":"FAIL");
    return ok?0:1;
}
