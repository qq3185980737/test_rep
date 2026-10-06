// test_deri0.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "deri0.h"
int main() {
    std::vector<double> e = {0,0,1,3,5};  // 1-based
    std::vector<double> diag(10,0), scalar(10,0);
    std::vector<int> nbo = {0,2,0,2};
    deri0(e, 4, scalar, diag, 0.0, nbo);
    bool ok = (std::abs(diag[1]-1.5)<1e-9 && std::abs(diag[2]-2.5)<1e-9 &&
               std::abs(diag[3]-1.0)<1e-9 && std::abs(diag[4]-2.0)<1e-9);
    std::printf("diag=(%g,%g,%g,%g) %s\n", diag[1],diag[2],diag[3],diag[4], ok?"PASS":"FAIL");
    return ok?0:1;
}
