// test_volume.cpp
#include <cstdio>
#include <cmath>
#include "volume.h"
int main() {
    double v1[3]={3,4,0};
    double a=volume(v1,1);
    double v2[6]={1,0,0,0,1,0};
    double ar=volume(v2,2);
    double v3[9]={1,0,0,0,1,0,0,0,1};
    double vol=volume(v3,3);
    bool ok=(std::fabs(a-5)<1e-9 && std::fabs(ar-1)<1e-9 && std::fabs(vol-1)<1e-9);
    std::printf("len=%g area=%g vol=%g %s\n",a,ar,vol,ok?"PASS":"FAIL");
    return ok?0:1;
}
