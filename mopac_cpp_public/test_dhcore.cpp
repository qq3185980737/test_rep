// test_dhcore.cpp
#include <cstdio>
#include "dhcore.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
extern "C" void h1elec_(int*, int*, double*, double*, double* di) { for (int i=0;i<81;++i) di[i]=0; }
extern "C" void rotate_(int*, int*, double*, double*, double*, int* kr, double*, double*, double* en) { *en=0; }
int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    numat = 2; norbs = 2;
    nat.resize(3); nfirst.resize(3); nlast.resize(3);
    nat[1]=1; nat[2]=1; nfirst[1]=1; nfirst[2]=2; nlast[1]=1; nlast[2]=2;
    double coord[6]={0,0,0, 1,0,0}, h[10], ww[10], e;
    dhcore(coord, h, ww, e, 1, 1, 0.01);
    std::printf("dhcore e=%g PASS\n", e);
    return 0;
}
