// test_fmat.cpp
#include <cstdio>
#include "fmat.h"
int main() {
    // f(x)=x^2, g=2x. delta=0.1, +0.5d=0.05 g=0.1; -0.5d=-0.05 g=-0.1.
    double h = fmat_second(0.1, -0.1, 0.1, 1.0);
    bool ok=(h==1.0);
    std::printf("fmat_second=%g %s\n",h, ok?"PASS":"FAIL");
    return ok?0:1;
}
