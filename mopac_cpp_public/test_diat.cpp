// test_diat.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "diat.h"
int main() {
    double s = ss_overlap(1,1,1,1,1, 1.0,1.0, 1.0, 1.0);
    bool ok = std::isfinite(s);
    std::printf("ss(1s,1s r=1)=%g %s\n", s, ok?"PASS":"FAIL");
    return ok?0:1;
}
