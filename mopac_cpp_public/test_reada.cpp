// test_reada.cpp
#include <cstdio>
#include <cmath>
#include "reada.h"
int main() {
    double v=reada(" abc1.23e2 x",4);
    bool ok=std::fabs(v-123)<1e-9;
    std::printf("v=%g %s\n",v,ok?"PASS":"FAIL");
    return ok?0:1;
}
