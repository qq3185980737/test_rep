// test_conref_C.cpp
#include <cstdio>
#include "conref_C.h"
int main() {
    bool ok = (conref_C::fpcref[1][3] == 0.5291772083 &&
               conref_C::fpcref[1][4] == 27.2113834 &&
               conref_C::fpcref[2][3] == 0.529167);
    std::printf("a0_98=%.10f eVau=%.7f a0_old=%.6f %s\n",
                conref_C::fpcref[1][3], conref_C::fpcref[1][4],
                conref_C::fpcref[2][3], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
