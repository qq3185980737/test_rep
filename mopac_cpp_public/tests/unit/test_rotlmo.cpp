// test_rotlmo.cpp
#include <cstdio>
#include <vector>
#include "rotlmo.h"
#include "MOZYME_C.h"
#include "molkst_C.h"
using namespace molkst_C;
using namespace MOZYME_C;
int main() {
    norbs=4; nelecs=2;
    iorbs.assign(3,4);
    nncf.assign(2,0); ncf.assign(2,1); ncocc.assign(2,0); icocc.assign(2,0);
    cocc.assign(8,0.0);
    nnce.assign(4,0); nce.assign(4,1); ncvir.assign(4,0); icvir.assign(4,0);
    cvir.assign(16,0.0);
    double rv[9]={1,0,0, 0,1,0, 0,0,1};
    rotlmo(rv);
    std::printf("rotlmo PASS\n");
    return 0;
}
