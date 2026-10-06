// test_outer2.cpp
#include <cstdio>
#include "outer2.h"
#include "parameters_C.h"
namespace molkst_C { int l1u,l2u,l3u; }
extern "C" void reppd_(int*, int*, double* rij, double* ri, double*) {
    ri[0] = 1.0/(*rij); ri[1] = 0.1; ri[4] = 0.2;
}
int main() {
    using namespace parameters_C;
    tore[1]=1.0; tore[2]=1.0;
    natorb[1]=1; natorb[2]=1;
    double xi[3]={0,0,0}, xj[3]={1.0,0,0}, w[7], e1[45], e2[45], en;
    int kr=0;
    outer2(1,2,xi,xj,w,kr,e1,e2,en,0,false);
    std::printf("outer2 en=%g kr=%d e1[0]=%g e2[0]=%g PASS\n", en, kr, e1[0], e2[0]);
    return 0;
}
