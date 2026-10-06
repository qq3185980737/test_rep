// test_mxm.cpp
#include <cstdio>
#include <cmath>
#include "mxm.h"
int main() {
    double a[4]={1,3,2,4},b[4]={5,7,6,8},c[4];
    mxm(a,2,b,2,c,2);
    bool ok=(std::fabs(c[0]-19)<1e-9 && std::fabs(c[1]-43)<1e-9 &&
             std::fabs(c[2]-22)<1e-9 && std::fabs(c[3]-50)<1e-9);
    std::printf("c=(%g,%g,%g,%g) %s\n",c[0],c[1],c[2],c[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
