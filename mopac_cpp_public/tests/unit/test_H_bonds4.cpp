// test_H_bonds4.cpp
#include <cstdio>
#include <cmath>
#include "H_bonds4.h"
int main() {
    double dp; double v=poly(1.0,true,dp); bool a=(std::fabs(v-25.4629)<0.01 && dp==0.0);
    double v2=poly(2.0,true,dp); bool b=(v2>0);
    std::printf("poly(1)=%g poly(2)=%g %s\n",v,v2,(a&&b)?"PASS":"FAIL");
    return (a&&b)?0:1;
}
