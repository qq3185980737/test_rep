// test_esp_C.cpp — verify esp_C module data (tf/fac/ixn..jzn) against Fortran.
#include "esp_C.h"

#include <cstdio>

using namespace esp_C;

int main() {
    bool ok = true;
    if (tf[0] != 33.0 || tf[1] != 37.0 || tf[2] != 41.0) { std::printf("FAIL tf\n"); ok = false; }
    double fac_expect[8] = {1, 1, 2, 6, 24, 120, 720, 5040};
    for (int i = 0; i < 8; ++i)
        if (fac[i] != fac_expect[i]) { std::printf("FAIL fac[%d]=%g\n", i, fac[i]); ok = false; }
    // Fortran: ixn = / 0, 4,0,0, 8,0,0, 4,4,0 /
    int ixn_expect[11] = {0, 0, 4, 0, 0, 8, 0, 0, 4, 4, 0};
    int iyn_expect[11] = {0, 0, 0, 4, 0, 0, 8, 0, 4, 0, 4};
    int izn_expect[11] = {0, 0, 0, 0, 4, 0, 0, 8, 0, 4, 4};
    int jxn_expect[11] = {0, 0, 1, 0, 0, 2, 0, 0, 1, 1, 0};
    int jyn_expect[11] = {0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 1};
    int jzn_expect[11] = {0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1};
    for (int i = 0; i < 11; ++i) {
        if (ixn[i] != ixn_expect[i] || iyn[i] != iyn_expect[i] || izn[i] != izn_expect[i] ||
            jxn[i] != jxn_expect[i] || jyn[i] != jyn_expect[i] || jzn[i] != jzn_expect[i]) {
            std::printf("FAIL index table[%d]: %d %d %d %d %d %d\n", i,
                        ixn[i], iyn[i], izn[i], jxn[i], jyn[i], jzn[i]);
            ok = false;
            break;
        }
    }
    // fv/dex start zeroed; tf entry 2 used for ESP fitting exponents
    if (fv[0][0] != 0.0 || dex[0] != 0.0) { std::printf("FAIL zero init\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
