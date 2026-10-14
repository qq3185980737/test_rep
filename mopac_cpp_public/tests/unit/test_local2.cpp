// test_local2.cpp
#include <cstdio>
#include <cmath>
#include "local2.h"
int main() {
    double c[4]={1,0,0,1}; int nf[2]={1,1},nl[2]={2,2};
    local2(c,2,2,nf,nl,2);
    bool ok=(std::fabs(c[0]-1)<1e-9 && std::fabs(c[3]-1)<1e-9);
    std::printf("local2 c0=%g c3=%g %s\n",c[0],c[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
