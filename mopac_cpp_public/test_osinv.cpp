// test_osinv.cpp
#include <cstdio>
#include <cmath>
#include "osinv.h"
int main() {
    double a[4]={2,0,0,3},d;
    osinv(a,2,d);
    bool ok=(std::fabs(a[0]-0.5)<1e-9 && std::fabs(a[3]-1.0/3)<1e-9);
    std::printf("osinv a0=%g a3=%g d=%g %s\n",a[0],a[3],d,ok?"PASS":"FAIL");
    return ok?0:1;
}
