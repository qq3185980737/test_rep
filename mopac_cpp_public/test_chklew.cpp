// test_chklew.cpp
#include <cstdio>
#include <vector>

#include "chklew.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"

int main() {
    using namespace MOZYME_C;
    Lewis_tot = 0;
    Lewis_elem.assign(3, std::vector<int>(10, 0));
    iz = {0, 4, 4};
    ib = {0, 4, 4};
    ions = {0, 0, 0};
    int t = 0;
    add_Lewis_element(1, 2, 0, t);   // sigma bond
    bool ok = (Lewis_tot == 1 && iz[1] == 3 && iz[2] == 3 && ib[1] == 3 && ib[2] == 3 && t == 1 &&
               Lewis_elem[1][1] == 1 && Lewis_elem[2][1] == 2);
    std::printf("sigma bond: t=%d iz=(%d,%d) ib=(%d,%d) %s\n", t, iz[1], iz[2], ib[1], ib[2], ok ? "PASS" : "FAIL");
    add_Lewis_element(1, 0, 0, t);   // lone pair on atom 1
    bool ok2 = (iz[1] == 1 && ib[1] == 2 && t == 2 && ions[1] == 0);
    std::printf("lone pair: t=%d iz1=%d ib1=%d %s\n", t, iz[1], ib[1], ok2 ? "PASS" : "FAIL");
    return (ok && ok2) ? 0 : 1;
}
