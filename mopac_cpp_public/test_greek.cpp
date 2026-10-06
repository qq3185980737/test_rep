// test_greek.cpp
#include <cstdio>
#include "greek.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
int main() {
    using namespace common_arrays_C;
    nat.assign(3,0); nbonds.assign(3,0);
    ibonds.assign(3, std::vector<int>(3,0));
    txtatm.assign(3, std::string(26,' '));
    MOZYME_C::nbackb.assign(4,0);
    greek(1);
    std::printf("greek PASS\n"); return 0;
}
