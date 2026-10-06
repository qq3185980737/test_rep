// test_check.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "check.h"
#include "MOZYME_C.h"

int main() {
    MOZYME_C::ws.assign(2, 0.0);
    std::vector<int> nnc = {0, 0}, nc = {0, 1}, icvec = {0, 1};
    std::vector<int> iorbs = {0, 1};
    std::vector<int> ncvec = {0, 0};
    std::vector<double> cvec = {0, 0.5};

    check(1, nnc, nc, icvec, iorbs, ncvec, cvec);
    bool ok = std::fabs(cvec[1] - 1.0) < 1e-9 && std::fabs(MOZYME_C::ws[1] - 0.25) < 1e-9;
    std::printf("cvec=%.3f ws=%.3f %s\n", cvec[1], MOZYME_C::ws[1], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
