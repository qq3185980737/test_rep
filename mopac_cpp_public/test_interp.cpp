// test_interp.cpp
#include <cstdio>
#include <cmath>
#include "interp.h"
int main() {
    // f(x)=(x-0.5)^2 sampled at 0 and 1
    double x[12]={0,1},f[12]={0.25,0.25},df[12]={-1,1};
    double xmin;
    spline(x,f,df,2.0,-2.0,xmin,2);
    bool ok=(std::fabs(xmin-0.5)<0.25);
    std::printf("spline xmin=%g %s\n",xmin,ok?"PASS":"FAIL");
    return ok?0:1;
}
