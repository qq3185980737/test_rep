// test_babbbc.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "babbbc.h"
#include "meci_C.h"

int main() {
    const int nmos = 2;
    std::vector<int> iocca1 = {0, 1, 0};   // alpha in MO1
    std::vector<int> ioccb1 = {0, 1, 0};  // state1: beta in MO1
    std::vector<int> ioccb2 = {0, 0, 1};  // state2: beta in MO2
    meci_C::occa = {0.0, 1.0, 1.0};

    // xy(p,q,r,s), 1-based. Only xy(1,2,2,2)=3.0 matters.
    std::vector<std::vector<std::vector<std::vector<double>>>> xy(
        nmos + 1, std::vector<std::vector<std::vector<double>>>(
                      nmos + 1, std::vector<std::vector<double>>(
                                    nmos + 1, std::vector<double>(nmos + 1, 0.0))));
    xy[1][2][2][2] = 3.0;

    double val = babbbc(iocca1, ioccb1, ioccb2, nmos, xy);
    // Expected: k=1 terms cancel; k=2 gives -xy(1,2,2,2) = -3.0; ij=0 even.
    bool ok = std::fabs(val - (-3.0)) < 1e-12;
    std::printf("babbbc=%.4f (expect -3.0) %s\n", val, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
