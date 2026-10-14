// test_dist2.cpp
#include <cstdio>
#include <array>
#include "dist2.h"
int main() {
    std::array<double,3> a={0,0,0}, b={1,2,2};
    double d=dist2(a,b), dt=dot1(a,b);
    bool ok=(d==9.0 && dt==0.0);
    std::printf("dist2=%g dot1=%g %s\n",d,dt,ok?"PASS":"FAIL");
    return ok?0:1;
}
