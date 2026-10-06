// test_deri22.cpp
#include <cstdio>
#include <vector>
#include "deri22.h"
int main() {
    std::vector<std::vector<double>> c, work;
    std::vector<double> foc2, ab(2,0), w, diag(2,1), scalar(2,3), fci;
    std::vector<double> b = {0,2.0};
    deri22(c, b, work, foc2, ab, 1, fci, w, diag, scalar, 1);
    bool ok = (b[1] == 2.0);
    std::printf("b after un/rescale=%g %s\n", b[1], ok?"PASS":"FAIL");
    return ok?0:1;
}
