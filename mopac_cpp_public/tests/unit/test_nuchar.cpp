// test_nuchar.cpp
#include <cstdio>
#include <cmath>
#include "nuchar.h"
int main() {
    char s[32]=" 1.23  4.5e2  ";
    double v[41]; int n;
    nuchar(s,13,v,n);
    bool ok=(n==2 && std::fabs(v[1]-1.23)<1e-9 && std::fabs(v[2]-450)<1e-9);
    std::printf("n=%d v1=%g v2=%g %s\n",n,v[1],v[2],ok?"PASS":"FAIL");
    return ok?0:1;
}
