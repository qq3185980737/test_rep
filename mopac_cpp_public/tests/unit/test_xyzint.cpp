// test_xyzint.cpp
#include <cstdio>
#include "xyzint.h"
#include "funcon_C.h"
namespace funcon_C { double pi = 3.14159265358979; }
extern "C" void bangle_(double*, int*, int*, int*, double* s) { *s = 0.1; }
extern "C" void dihed_(double*, int*, int*, int*, int*, double* d) { *d = 0.0; }
int main() {
    double xyz[12], geo[12];
    int na[5]={0,0,0,0,0}, nb[5]={0,0,0,0,0}, nc[5]={0,0,0,0,0};
    // 4 atoms
    for (int i=0;i<12;++i) xyz[i]=0;
    xyz[0]=0; xyz[3]=1; xyz[6]=2; xyz[9]=3;
    xyz[1]=0; xyz[4]=0; xyz[7]=0; xyz[10]=0;
    xyz[2]=0; xyz[5]=0; xyz[8]=0; xyz[11]=0;
    xyzint(xyz, 4, na, nb, nc, 57.29578, geo);
    std::printf("xyzint r1=%g na2=%d PASS\n", geo[0], na[2]);
    return 0;
}
