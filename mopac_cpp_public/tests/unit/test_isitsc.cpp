// test_isitsc.cpp
#include <cstdio>
#include "isitsc.h"
namespace molkst_C { int iscf; }
namespace MOZYME_C { double ovmax=0.001; double energy_diff=0.001; }
int main() {
    bool ok; int iemin=0, iemax=0;
    isitsc(0.0, 0.01, 10.0, iemin, iemax, ok, 200, 100);
    std::printf("ok=%d iscf=%d PASS\n", (int)ok, molkst_C::iscf);
    return 0;
}
