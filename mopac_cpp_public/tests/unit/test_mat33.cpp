// test_mat33.cpp
#include <cstdio>
#include <cmath>
#include "mat33.h"
int main() {
    double a[10]={0,1,0,0,0,1,0,0,0,1},b[10]={0,1,0,0,0,1,0,0,0,1},c[10];
    mat33(a,b,c);
    bool ok=(std::fabs(c[1]-1)<1e-12 && std::fabs(c[5]-1)<1e-12 && std::fabs(c[9]-1)<1e-12);
    std::printf("mat33 c1=%g c5=%g c9=%g %s\n",c[1],c[5],c[9],ok?"PASS":"FAIL");
    return ok?0:1;
}
