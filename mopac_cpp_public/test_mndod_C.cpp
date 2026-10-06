// test_mndod_C.cpp — verify mndod_C tables: iii/iiid/intij/intkl/intrep counts,
// spot values, and empty runtime arrays usable by mndod.cpp.
#include "mndod_C.h"

#include <cstdio>

using namespace mndod_C;

int main() {
    bool ok = true;
    if ((int)iii.size() != 108 || (int)iiid.size() != 108) { std::printf("FAIL iii/iiid size\n"); ok = false; }
    if ((int)intij.size() != 244 || (int)intkl.size() != 244 || (int)intrep.size() != 244) {
        std::printf("FAIL 243-tables size: %d %d %d\n", (int)intij.size(), (int)intkl.size(), (int)intrep.size());
        ok = false;
    }
    if ((int)isym.size() != 492) { std::printf("FAIL isym size\n"); ok = false; }
    // Fortran: iii = 2*1, 8*2, 8*3, 18*4, 18*5, 32*6, 21*0
    if (iii[1] != 1 || iii[2] != 1 || iii[3] != 2 || iii[10] != 2 || iii[11] != 3 ||
        iii[19] != 4 || iii[37] != 5 || iii[55] != 6 || iii[107] != 0) {
        std::printf("FAIL iii values\n"); ok = false;
    }
    // Fortran: iiid = 30*3, 18*4, 32*5, 6*6, 21*0
    if (iiid[1] != 3 || iiid[30] != 3 || iiid[31] != 4 || iiid[48] != 4 ||
        iiid[49] != 5 || iiid[80] != 5 || iiid[81] != 6 || iiid[86] != 6 || iiid[87] != 0) {
        std::printf("FAIL iiid values\n"); ok = false;
    }
    // intij/intkl/intrep first and last entries (Fortran data order)
    if (intij[1] != 1 || intkl[1] != 15 || intrep[1] != 1) { std::printf("FAIL table[1]\n"); ok = false; }
    if (intij[243] != 45 || intkl[243] != 45 || intrep[243] != 29) {
        std::printf("FAIL table[243]: %d %d %d\n", intij[243], intkl[243], intrep[243]);
        ok = false;
    }
    // runtime arrays are preallocated and writable
    indexd.resize(10, std::vector<int>(10, 0));
    indexd[5][5] = 7;
    if (indexd[5][5] != 7) { std::printf("FAIL indexd\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
