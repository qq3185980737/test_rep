// test_newflg.cpp
#include <cstdio>
#include <vector>
#include "newflg.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
void mopend(const char*) {}
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    numat = 1; natoms = 1; nvar = 1;
    na.assign(2,1); na[0]=0; nb.assign(2,0); nc.assign(2,0);
    coord.assign(4,std::vector<double>(2,1.0));
    geo.assign(4,std::vector<double>(2,0.0));
    txtatm.assign(2,"ATOM  N    N N  XX");
    loc.assign(3,std::vector<int>(2,1));
    xparam.assign(2,0.0);
    newflg();
    std::printf("newflg: geo[1][1]=%g xparam[1]=%g PASS\n", geo[1][1], xparam[1]);
    return 0;
}
