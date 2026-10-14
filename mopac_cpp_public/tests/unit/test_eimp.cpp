// test_eimp.cpp
#include <cstdio>
#include "eimp.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numat=1;
    eimp();
    std::printf("eimp PASS\n");
    return 0;
}
