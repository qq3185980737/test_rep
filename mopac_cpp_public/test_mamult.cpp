// test_mamult.cpp
#include <cstdio>
#include <cmath>
#include "mamult.h"
int main() {
    double a[4]={0,1,0,2},b[4]={0,1,0,2},c[4]={0,0,0,0};
    mamult(a,b,c,2,0.0);
    // c1=a11*b11+a21*b21=1; c2=a21*b11+a22*b21=0; c3=a21*b21+a22*b22=4
    bool ok=(std::fabs(c[1]-1)<1e-12 && std::fabs(c[2])<1e-12 && std::fabs(c[3]-4)<1e-12);
    std::printf("c1=%g c2=%g c3=%g %s\n",c[1],c[2],c[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
