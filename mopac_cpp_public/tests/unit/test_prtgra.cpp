// test_prtgra.cpp
#include <cstdio>
#include "prtgra.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int iw=6;
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    natoms=1; nvar=3; maxtxt=0; na1=0;
    grad.assign(4,0.5); xparam.assign(4,0.1);
    labels.assign(2,1); na.assign(2,0);
    loc.assign(3,std::vector<int>(4,0));
    loc[1][1]=1; loc[2][1]=1;
    loc[1][2]=1; loc[2][2]=2;
    loc[1][3]=1; loc[2][3]=3;
    prtgra();
    std::printf("prtgra PASS\n");
    return 0;
}
