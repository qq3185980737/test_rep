// test_molmec_C.cpp — molmec_C module: nhco(4,4000), nnhco, htype.
#include "molmec_C.h"

#include <cstdio>

using namespace molmec_C;

int main() {
    bool ok = true;
    if (nnhco != 0 || htype != 0.0) { std::printf("FAIL defaults\n"); ok = false; }
    nnhco = 3;
    nhco[2][3] = 42; nhco[1][4000] = 7;
    if (nnhco != 3 || nhco[2][3] != 42 || nhco[1][4000] != 7) { std::printf("FAIL writes\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
