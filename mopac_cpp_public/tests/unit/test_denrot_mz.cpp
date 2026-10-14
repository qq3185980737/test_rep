// test_denrot_mz.cpp
#include <cstdio>
#include "denrot_for_MOZYME.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numat = 0;
    denrot_for_MOZYME();
    std::printf("denrot_for_MOZYME(empty) PASS\n");
    return 0;
}
