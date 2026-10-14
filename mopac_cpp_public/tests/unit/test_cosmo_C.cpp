// test_cosmo_C.cpp
#include <cstdio>
#include "cosmo_C.h"
int main() {
    bool ok = (cosmo_C::nppa == 1082 && !cosmo_C::iseps && cosmo_C::solv_energy == 0.0);
    std::printf("nppa=%d iseps=%d solv=%g %s\n",
                cosmo_C::nppa, (int)cosmo_C::iseps, cosmo_C::solv_energy, ok?"PASS":"FAIL");
    return ok?0:1;
}
