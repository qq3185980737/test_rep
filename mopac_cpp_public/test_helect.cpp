// test_helect.cpp
#include <cstdio>
#include <cmath>
#include "helect.h"
int main() {
    double p[4]={0,1,2,3},h[4]={0,1,0,1},f[4]={0,0,0,0};
    double v=helect(2,p,h,f);
    double expect=2*0+0.5*(1*1+3*1);
    bool ok=(std::fabs(v-expect)<1e-12);
    std::printf("helect=%g expect=%g %s\n",v,expect,ok?"PASS":"FAIL");
    return ok?0:1;
}
