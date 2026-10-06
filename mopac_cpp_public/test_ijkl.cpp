// test_ijkl.cpp
#include <cstdio>
#include "ijkl.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
namespace cosmo_C { bool iseps=false; }
int main() {
    using namespace common_arrays_C;
    molkst_C::norbs=4; molkst_C::numat=2;
    nfirst.assign(3,0); nlast.assign(3,0);
    nfirst[1]=1; nlast[1]=2; nfirst[2]=3; nlast[2]=4;
    w.assign(20,1.0);
    double cp[8]={1,0,0,0, 0,1,0,0};
    double cf[16]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    double dijkl[4*2*3];
    double cij[20],ckl[20],wcij[20];
    double xy[2*2*2*2];
    ijkl(cp,cf,2,2,dijkl,cij,ckl,wcij,xy);
    std::printf("ijkl PASS\n");
    return 0;
}
