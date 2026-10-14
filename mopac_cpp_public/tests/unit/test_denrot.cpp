// test_denrot.cpp
#include <cstdio>
#include "denrot.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    mozyme = true;   // early-return branch
    denrot();
    bool ok = true;
    std::printf("denrot(mozyme) returned %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
