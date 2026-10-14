// test_H_bond_correction_bits.cpp
#include <cstdio>
#include <cmath>
#include "H_bond_correction_bits.h"
int main() {
    double lo=truncation(0.0,3.0,0.5), mid=truncation(3.0,3.0,0.5), hi=truncation(10.0,3.0,0.5);
    bool ok=(std::fabs(lo-3.0)<1e-9 && std::fabs(mid-3.125)<1e-9 && std::fabs(hi-10.0)<1e-9);
    std::printf("trunc lo=%g mid=%g hi=%g %s\n",lo,mid,hi,ok?"PASS":"FAIL");
    return ok?0:1;
}
