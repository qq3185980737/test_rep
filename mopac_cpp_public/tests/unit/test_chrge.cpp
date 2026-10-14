// test_chrge.cpp
#include <cstdio>
#include <vector>
#include "chrge.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
namespace common_arrays_C {
    std::vector<int> nfirst, nlast;
}
namespace molkst_C {
    int numat, norbs; bool mozyme;
}
using namespace common_arrays_C;
using namespace molkst_C;
int main() {
    numat = 2; norbs = 3; mozyme = false;
    nfirst = {0, 1, 3}; nlast = {0, 2, 3};
    std::vector<double> pin(7, 0.0);
    pin[1] = 1.0; pin[3] = 1.0; pin[6] = 1.0;
    std::vector<double> qout(3, -1.0);
    chrge(pin, qout);
    bool ok = (qout[1] == 2.0 && qout[2] == 1.0);
    std::printf("q=(%g,%g) %s\n", qout[1], qout[2], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
