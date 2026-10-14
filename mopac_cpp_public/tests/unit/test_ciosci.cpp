// test_ciosci.cpp
#include <cstdio>
#include <vector>

#include "ciosci.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

int main() {
    using namespace common_arrays_C;
    using namespace funcon_C;
    using namespace meci_C;
    using namespace molkst_C;
    using namespace parameters_C;
    numat = 1; norbs = 1; nmos = 1; nstate = 1; lab = 1;
    nfirst = {0, 1}; nlast = {0, 1}; nat = {0, 1};
    coord.assign(4, std::vector<double>(2, 0.0));
    dd[1] = 0.0; a0 = 0.529;
    microa.assign(2, std::vector<int>(2, 0));
    microb.assign(2, std::vector<int>(2, 0));
    std::vector<std::vector<double>> v(2, std::vector<double>(2, 0.0));
    v[1][1] = 1.0;
    std::vector<double> cin(2, 0.0); cin[1] = 1.0;
    std::vector<std::vector<double>> osc;
    ciosci(v, 1, osc, cin);
    bool ok = (osc[1][1] == 0.0 && osc[2][1] == 0.0 && osc[3][1] == 0.0);
    std::printf("oscil=(%g,%g,%g) %s\n", osc[1][1], osc[2][1], osc[3][1], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
