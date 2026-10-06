// test_wrttxt.cpp
#include <cstdio>
#include "wrttxt.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    refkey.assign(6, " NULL");
    refkey[0] = "TEST";
    wrttxt(6);
    std::printf("wrttxt links OK PASS\n");
    return 0;
}
