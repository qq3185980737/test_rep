// test_schmit.cpp
#include <cstdio>
#include <cmath>
#include "schmit.h"
int main() {
    double u[9]={1,0,0, 1,1,0, 0,1,0};
    schmit(u,3,3);
    auto col=[&](int c){ return &u[(c-1)*3]; };
    auto dot=[&](int a,int b){ double s=0; for(int i=0;i<3;++i) s+=col(a)[i]*col(b)[i]; return s; };
    double n1=dot(1,1),n2=dot(2,2),n3=dot(3,3),d12=dot(1,2),d13=dot(1,3),d23=dot(2,3);
    bool ok=(std::fabs(n1-1)<1e-9&&std::fabs(n2-1)<1e-9&&std::fabs(n3-1)<1e-9&&
             std::fabs(d12)<1e-9&&std::fabs(d13)<1e-9&&std::fabs(d23)<1e-9);
    std::printf("n=(%g,%g,%g) d=(%g,%g,%g) %s\n",n1,n2,n3,d12,d13,d23,ok?"PASS":"FAIL");
    return ok?0:1;
}
