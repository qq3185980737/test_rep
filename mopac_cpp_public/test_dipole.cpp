// test_dipole.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dipole.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    numat=1; numcal=1; keywrd="";
    nfirst.assign(2,1); nlast.assign(2,1); nat.assign(2,1);
    q.assign(2,0.0); q[1]=1.0;
    std::vector<double> pv(10,0.0), dv(3);
    std::vector<std::vector<double>> crd(3, std::vector<double>(2,0.0));
    crd[0][1]=1.0;
    double d = dipole(pv, crd, dv, 0);
    bool ok = std::isfinite(d);
    std::printf("dipole=%g %s\n", d, ok?"PASS":"FAIL");
    return ok?0:1;
}
