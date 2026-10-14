// test_dfock2.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dfock2.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    norbs=2; numat=2;
    nfirst = {0,1,1}; nlast = {0,1,1};
    std::vector<double> fv(3,0), ptot={0,10.0}, pv={0,4.0}, wv={0,2.0};
    dfock2(fv, ptot, pv, wv, 2, 1);
    double expect = 2.0*10.0*2.0 - 4.0*2.0;
    bool ok = (std::abs(fv[1]-expect)<1e-9);
    std::printf("fv[1]=%g expect=%g %s\n", fv[1], expect, ok?"PASS":"FAIL");
    return ok?0:1;
}
