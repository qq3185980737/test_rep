// test_h1elec.cpp
#include <cstdio>
#include <cstring>
#include "h1elec.h"
#include "parameters_C.h"
#include "MOZYME_C.h"
namespace MOZYME_C { double cutofs = 100.0; }
extern "C" void diat_(int*, int*, double*, double* s) {
    for (int i=0;i<81;++i) s[i]=0.0;
    s[0] = 1.0;
}
int main() {
    using namespace parameters_C;
    for (int i=0;i<=N_PARAM;++i) { betas[i]=1.0; betap[i]=1.0; betad[i]=0.0; natorb[i]=1; }
    natorb[6] = 4; natorb[1] = 1;
    double xi[3]={0,0,0}, xj[3]={1.0,0,0}, s[81];
    h1elec(6,1,xi,xj,s);
    std::printf("h1elec s[0]=%g PASS\n", s[0]);
    return 0;
}
