// test_compfg.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "compfg.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "parameters_C.h"

int main() {
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace parameters_C;
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][2] = 1.0;      // r in bohr along x
    tore[1] = 1.0;
    ev = 27.2114; a0 = 0.529177; fpc_9 = 23.0605;
    po.assign(11, std::vector<double>(110, 0.0));
    // abond[1][1] left 0 => generic branch
    double e = xfac_value();
    bool ok = (e > 0.0 && std::isfinite(e));
    std::printf("xfac_value(H-H, r=1bohr)=%.3f %s\n", e, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
