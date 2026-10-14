// test_mkl_bits.cpp
#include <cstdio>
#include <cmath>
#include "mkl_bits.h"
int main() {
    double x[3]={1,2,3},y[3]={4,5,6};
    double d=ddot(3,x,1,y,1);
    bool ok=(std::fabs(d-32)<1e-9);
    std::printf("ddot=%g %s\n",d,ok?"PASS":"FAIL");
    return ok?0:1;
}
