// test_reorth.cpp
#include <cstdio>
#include "reorth.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    molkst_C::numat=1; molkst_C::norbs=1;
    using namespace MOZYME_C;
    nvirtual=0; noccupied=0;
    common_arrays_C::nfirst.assign(2,1); iorbs.assign(2,1);
    double wsv[5]; reorth(wsv);
    std::printf("reorth PASS\n"); return 0;
}
