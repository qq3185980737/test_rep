// test_dipole_for_MOZYME.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dipole_for_MOZYME.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
int main() {
    using namespace molkst_C; using namespace common_arrays_C;
    numat=1; numcal=1; keywrd="";
    nat.assign(2,1); q.assign(2,0.0); q[1]=1.0;
    parameters_C::ams.assign(108,1.0);
    p.assign(50,0.0);
    MOZYME_C::iorbs.assign(2,1);
    coord.assign(3, std::vector<double>(2,0.0));
    coord[0][1]=1.0;
    std::vector<double> dv(3);
    double d = dipole_for_MOZYME(dv,0);
    bool ok = std::isfinite(d);
    std::printf("dipole_MOZYME=%g %s\n", d, ok?"PASS":"FAIL");
    return ok?0:1;
}
