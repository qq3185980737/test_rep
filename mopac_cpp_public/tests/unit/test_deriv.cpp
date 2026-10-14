// test_deriv.cpp
#include <cstdio>
#include <vector>
#include "deriv.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    nvar = 1; numcal = 1; keywrd = "";
    std::vector<std::vector<double>> geo(4, std::vector<double>(4,0));
    std::vector<double> gradnt = {0, 2.0};
    deriv(geo, gradnt);
    bool ok = (gradnt[1] == 2.0 * 1e7);
    std::printf("gradnt=%g %s\n", gradnt[1], ok?"PASS":"FAIL");
    return ok?0:1;
}
