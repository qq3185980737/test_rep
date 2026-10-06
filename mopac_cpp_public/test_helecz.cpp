// test_helecz.cpp
#include <cstdio>
#include <vector>
#include "helecz.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
namespace molkst_C { int numat; }
namespace MOZYME_C { std::vector<int> iorbs; }
namespace common_arrays_C {
    std::vector<double> p, h, f;
}
extern "C" int ijbo_(int* i, int* j) { return (*i==1 && *j==1) ? 0 : -1; }
int main() {
    using namespace molkst_C;
    using namespace MOZYME_C;
    using namespace common_arrays_C;
    numat = 1;
    iorbs = {0, 1};
    p = {0, 2.0}; h = {0, -1.0}; f = {0, -1.0};
    double v = helecz();
    // diagonal: ed += p(1)*(h(1)+f(1)) = 2*(-2) = -4; ee = 0.5*ed = -2
    bool ok = (v == -2.0);
    std::printf("helecz=%g %s\n", v, ok?"PASS":"FAIL");
    return ok?0:1;
}
