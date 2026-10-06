// test_delsta.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "delsta.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
int main() {
    using namespace funcon_C;
    using namespace molkst_C;
    using namespace parameters_C;
    cutofp = 10.0; ev = 27.2114; fpc_9 = 23.0605; tore[1] = 1.0;
    std::vector<int> nat = {0,1,1}, iorbs = {0,1,1};
    std::vector<double> p(3, 0.0);
    std::vector<std::vector<double>> cdi(4, std::vector<double>(3, 0.0));
    cdi[1][1] = 1.0;  // atom1 at x=1, atom2 at 0
    std::vector<double> dstat(4, 0.0);
    delsta(nat, iorbs, p, cdi, dstat, 1, 2);
    double expect = -0.5 * fpc_9 * ev;  // q=1, rij=1, vect=1
    bool ok = (std::abs(dstat[1] - expect) < 1e-6 && std::abs(dstat[2]) < 1e-9);
    std::printf("dstat=(%g,%g,%g) expect=%g %s\n",
                dstat[1], dstat[2], dstat[3], expect, ok?"PASS":"FAIL");
    return ok?0:1;
}
