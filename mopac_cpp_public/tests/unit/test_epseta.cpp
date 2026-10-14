// test_epseta.cpp
#include <cstdio>
#include "epseta.h"
int main() {
    double e, a; epseta(e,a);
    bool ok=(e>0 && a>0);
    std::printf("epseta eps=%g eta=%g %s\n",e,a, ok?"PASS":"FAIL");
    return ok?0:1;
}
