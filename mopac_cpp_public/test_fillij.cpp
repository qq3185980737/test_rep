// test_fillij.cpp
#include <cstdio>
#include "overlaps_C.h"
#include "fillij.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
namespace overlaps_C { double cutof1=0, cutof2=0; }
namespace common_arrays_C {
    std::vector<std::vector<double>> coord;
    std::vector<std::vector<double>> tvec;
}
namespace MOZYME_C {
    std::vector<int> iij, numij, ijall, iijj;
}
extern "C" double reada_(const char*, int*, int) { return 1.0; }
int main() {
    using namespace molkst_C;
    numat = 50;
    natoms = 50;
    cutofp = 0;
    MOZYME_C::iorbs.assign(numat + 1, 1);
    common_arrays_C::coord.resize(3);
    for (int d = 0; d < 3; ++d) common_arrays_C::coord[d].assign(numat + 1, 0.0);
    fillij(true);
    bool ok = (overlaps_C::cutof1 >= overlaps_C::cutof2);
    std::printf("fillij cutof1=%g cutof2=%g mpack=%d n2elec=%d %s\n",
        overlaps_C::cutof1, overlaps_C::cutof2, mpack, n2elec, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
