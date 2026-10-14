// test_mtxm.cpp
#include <cstdio>
#include <cmath>
#include "mtxm.h"
int main() {
    double a[4]={1,3,2,4}; // a(2,2): col1=(1,3), col2=(2,4)
    double b[4]={5,7,6,8};
    double c[4];
    mtxm(a,2,b,2,c,2);
    bool ok=(std::fabs(c[0]-26)<1e-9 && std::fabs(c[1]-38)<1e-9 &&
             std::fabs(c[2]-30)<1e-9 && std::fabs(c[3]-44)<1e-9);
    std::printf("c=(%g,%g,%g,%g) %s\n",c[0],c[1],c[2],c[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
