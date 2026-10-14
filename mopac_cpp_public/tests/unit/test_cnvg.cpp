// test_cnvg.cpp
#include <cstdio>
#include <vector>

#include "cnvg.h"
#include "molkst_C.h"

int main() {
    using namespace molkst_C;
    norbs = 1; numcal = 1; keywrd = "";
    std::vector<double> pnew(2, 0.0), p(2, 0.0), p1(2, 0.0);
    pnew[1] = 2.0;
    double pl = -1.0;
    cnvg(pnew, p, p1, 1, pl);
    bool ok = (pl == 2.0 && pnew[1] == 2.0);
    std::printf("pl=%g pnew[1]=%g %s\n", pl, pnew[1], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
