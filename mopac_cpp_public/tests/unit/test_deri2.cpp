// test_deri2.cpp
#include <cstdio>
#include <vector>
#include "deri2.h"
int main() {
    std::vector<std::vector<double>> f, fd, fci;
    std::vector<double> dxyzr, diag, scalar, work;
    deri2(0, f, fd, fci, 0, 0, dxyzr, 1e-6, diag, scalar, work);
    std::printf("deri2(skeleton) PASS\n");
    return 0;
}
