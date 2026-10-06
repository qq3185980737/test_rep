// test_mxv.cpp
#include <cstdio>
#include <cmath>
#include "mxv.h"
int main() {
    double a[4]={1,3,2,4},x[3]={0,1,2},y[3];
    mxv(a,2,x,2,y);
    bool ok=(std::fabs(y[1]-5)<1e-9 && std::fabs(y[2]-11)<1e-9);
    std::printf("y=(%g,%g) %s\n",y[1],y[2],ok?"PASS":"FAIL");
    return ok?0:1;
}
