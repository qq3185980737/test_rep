// test_prtlmo.cpp
#include <cstdio>
#include "prtlmo.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
int iw = 6;

int main() {
    using namespace molkst_C; using namespace MOZYME_C; using namespace common_arrays_C;
    norbs = 2; noccupied = 1; nvirtual = 1;
    eigs.assign(3, 0.0);
    nat.assign(3, 1);
    iorbs.assign(3, 1);
    nncf.assign(3, 0); icocc.assign(3, 0); ncocc.assign(3, 0); cocc.assign(3, 0.0);
    ncf.assign(3, 0); nnce.assign(3, 0); icvir.assign(3, 0); ncvir.assign(3, 0);
    cvir.assign(3, 0.0); nce.assign(3, 0);
    prtlmo();
    std::printf("prtlmo PASS\n");
    return 0;
}
