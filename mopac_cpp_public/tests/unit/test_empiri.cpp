// test_empiri.cpp — M03 batch C: empiri (empirical formula) verification
// Expected strings follow real MOPAC output format:
//   "           Empirical Formula: C4 H22 N10 O2 Cu Cl6  =    45 atoms"
// count==1 elements print symbol only (no digit); count>1 digits are glued.
#include <cstdio>
#include <string>
#include "empiri.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

static int passed = 0, failed = 0;
static void check(const char* name, bool ok) {
    if (ok) { ++passed; std::printf("  [PASS] %s\n", name); }
    else { ++failed; std::printf("  [FAIL] %s\n", name); }
}

int main() {
    std::printf("M03 batch C tests (empiri)\n");
    // nat 1-based: nat[1..numat]
    {
        // C2 H3 O1 : single-count O -> no digit
        common_arrays_C::nat.assign(7, 0);
        common_arrays_C::nat[1] = 6; common_arrays_C::nat[2] = 6;
        common_arrays_C::nat[3] = 1; common_arrays_C::nat[4] = 1; common_arrays_C::nat[5] = 1;
        common_arrays_C::nat[6] = 8;
        molkst_C::numat = 6;
        empiri();
        check("empiri C2H3O single-count blank",
              molkst_C::formula == "           Empirical Formula: C2 H3 O  =     6 atoms");
    }
    {
        // Cu Cl2 : no C/H/N/O; encounter order Cl, Cu; Cl2 has digit, Cu single -> blank
        common_arrays_C::nat.assign(4, 0);
        common_arrays_C::nat[1] = 17; common_arrays_C::nat[2] = 17; common_arrays_C::nat[3] = 29;
        molkst_C::numat = 3;
        empiri();
        check("empiri CuCl2 two-letter + single",
              molkst_C::formula == "           Empirical Formula: Cl2 Cu  =     3 atoms");
    }
    {
        // CH4 : C1 single -> blank, H4 digit
        common_arrays_C::nat.assign(6, 0);
        common_arrays_C::nat[1] = 6;
        common_arrays_C::nat[2] = 1; common_arrays_C::nat[3] = 1; common_arrays_C::nat[4] = 1; common_arrays_C::nat[5] = 1;
        molkst_C::numat = 5;
        empiri();
        check("empiri CH4 leading C single",
              molkst_C::formula == "           Empirical Formula: C H4  =     5 atoms");
    }
    std::printf("M03 batch C: %d passed, %d failed\n", passed, failed);
    return failed ? 1 : 0;
}
