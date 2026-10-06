// test_mult.cpp
#include <cstdio>
#include <cmath>
#include "mult.h"
int main() {
    double c[4]={1,0,0,1},s[4]={1,0,0,1},v[4];
    mult(c,s,v,2);
    bool ok=(std::fabs(v[0]-1)<1e-9 && std::fabs(v[1])<1e-9 &&
             std::fabs(v[2])<1e-9 && std::fabs(v[3]-1)<1e-9);
    std::printf("v=(%g,%g,%g,%g) %s\n",v[0],v[1],v[2],v[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
