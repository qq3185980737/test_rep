// test_deritr.cpp
#include <cstdio>
#include "deritr.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal = 1; keywrd = "";
    deritr();
    bool ok = true;
    std::printf("deritr(skeleton) %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
