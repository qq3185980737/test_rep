// test_dfield.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dfield.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
namespace common_arrays_C {
    std::vector<int> nat;
    std::vector<double> p, dxyz;
}
namespace funcon_C { double ev, a0, fpc_9; }
namespace molkst_C { int numat; double efield[4]; }
namespace parameters_C { double tore[108]; }
int main() {
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace molkst_C;
    using namespace parameters_C;
    numat=1; ev=27.2114; a0=0.529177; fpc_9=23.0605;
    tore[1]=1.0; nat.resize(2,1); nat[1]=1;
    efield[1]=0; efield[2]=1.0; efield[3]=0;
    dxyz.assign(5, 0.0);
    p.assign(5,0.0);
    dfield();
    double expect = ev/a0*fpc_9;
    bool ok = (std::abs(dxyz[2]-expect)<1e-6 && std::abs(dxyz[1])<1e-9);
    std::printf("dxyz2=%g expect=%g %s\n", dxyz[2], expect, ok?"PASS":"FAIL");
    return ok?0:1;
}
