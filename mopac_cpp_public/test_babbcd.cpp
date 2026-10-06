// test_babbcd.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "babbcd.h"

int main() {
    const int nmos = 4;
    std::vector<int> iocca1 = {0, 0, 0, 0, 0};
    std::vector<int> ioccb1 = {0, 1, 1, 0, 0};  // state1 beta in (1,2)
    std::vector<int> iocca2 = {0, 0, 0, 0, 0};
    std::vector<int> ioccb2 = {0, 0, 0, 1, 1};  // state2 beta in (3,4)

    std::vector<std::vector<std::vector<std::vector<double>>>> xy(
        nmos + 1, std::vector<std::vector<std::vector<double>>>(
                      nmos + 1, std::vector<std::vector<double>>(
                                    nmos + 1, std::vector<double>(nmos + 1, 0.0))));
    xy[3][1][4][2] = 5.0;
    xy[3][2][4][1] = 2.0;

    double val = babbcd(iocca1, ioccb1, iocca2, ioccb2, nmos, xy);
    // i=3, j=4, k=1, l=2, ij=0 -> one=+1; result = xy(3,1,4,2)-xy(3,2,4,1) = 3.
    bool ok = std::fabs(val - 3.0) < 1e-12;
    std::printf("babbcd=%.4f (expect 3.0) %s\n", val, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
