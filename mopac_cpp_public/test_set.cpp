// test_set.cpp
#include <cstdio>
#include <cmath>
#include "set.h"
int main() {
    set_fn(1.0,1.0,1,2,2.0,1);
    double alpha=0.5*2.0*(1+1); // =2
    double expect_a1=std::exp(-alpha)/alpha;
    bool ok=std::fabs(overlaps_C::a[1]-expect_a1)<1e-9;
    std::printf("a1=%g expect=%g %s\n",overlaps_C::a[1],expect_a1,ok?"PASS":"FAIL");
    return ok?0:1;
}
