// test_prttim.cpp
#include <cstdio>
#include <cmath>
#include "prttim.h"
int main() {
    double t; char c;
    prttim(7200,t,c);
    bool ok=(c=='H' && std::fabs(t-2)<1e-9);
    std::printf("t=%g c=%c %s\n",t,c,ok?"PASS":"FAIL");
    return ok?0:1;
}
