// test_chkion.cpp
#include <cstdio>
#include <vector>

#include "chkion.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    numat = 1; nelecs = 0;
    nat = {0, 8};
    nbonds = {0, 0};
    ibonds.assign(4, std::vector<int>(4, 0));
    MOZYME_C::ions = {0, 0};
    MOZYME_C::ib = {0, 0};
    MOZYME_C::iz = {0, 0};
    MOZYME_C::Lewis_tot = 0;
    MOZYME_C::Lewis_elem.assign(3, std::vector<int>(1, 0));
    std::vector<char> ac(2, ' ');
    std::vector<int> ox(2, -9);
    int nlp = 0;

    chkion(ox, nlp, ac);
    bool ok = (ox[1] == 0);
    std::printf("ox_calc(O isolated)=%d %s\n", ox[1], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
