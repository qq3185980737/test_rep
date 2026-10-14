// test_deri1.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "deri1.h"
int main() {
    double grad = 0;
    std::vector<double> fv={0,2,4}, fd(5,0), scalar={0,3,3};
    std::vector<std::vector<double>> work(3, std::vector<double>(3,0));
    deri1(1, grad, fv, 2, fd, scalar, work);
    bool ok = (std::abs(fv[1]-6.0)<1e-9 && std::abs(fv[2]-12.0)<1e-9);
    std::printf("f=(%g,%g) %s\n", fv[1], fv[2], ok?"PASS":"FAIL");
    return ok?0:1;
}
