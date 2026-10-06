// test_minv.cpp
#include <cstdio>
#include <cmath>
#include "minv.h"
int main() {
    double a[4]={2,0,0,3}, d;
    minv(a,2,d);
    bool ok=(std::fabs(a[0]-0.5)<1e-9 && std::fabs(a[3]-1.0/3)<1e-9 && std::fabs(d-6)<1e-9);
    std::printf("a0=%g a3=%g det=%g %s\n",a[0],a[3],d,ok?"PASS":"FAIL");
    return ok?0:1;
}
