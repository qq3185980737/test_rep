// test_supdot.cpp
#include <cstdio>
#include <cmath>
#include "supdot.h"
int main() {
    double s[3], h[4]={0,1,2,3}, g[3]={0,1,1};
    supdot(s,h,g,2);
    bool ok=(std::fabs(s[1]-3)<1e-9 && std::fabs(s[2]-5)<1e-9);
    std::printf("s1=%g s2=%g %s\n",s[1],s[2],ok?"PASS":"FAIL");
    return ok?0:1;
}
