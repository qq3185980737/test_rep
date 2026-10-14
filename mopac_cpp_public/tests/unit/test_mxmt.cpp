// test_mxmt.cpp
#include <cstdio>
#include <cmath>
#include "mxmt.h"
int main() {
    double a[4]={1,3,2,4},b[4]={5,7,6,8},c[4];
    mxmt(a,2,b,2,c,2);
    bool ok=(std::fabs(c[0]-17)<1e-9 && std::fabs(c[1]-39)<1e-9 &&
             std::fabs(c[2]-23)<1e-9 && std::fabs(c[3]-53)<1e-9);
    std::printf("c=(%g,%g,%g,%g) %s\n",c[0],c[1],c[2],c[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
