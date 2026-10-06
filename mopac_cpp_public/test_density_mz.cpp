// test_density_mz.cpp
#include <cstdio>
#include <vector>
#include "density_for_MOZYME.h"
int main() {
    std::vector<double> p(5, 9.0), partp(5, 1.0);
    density_for_MOZYME(p, 0, 0, partp);
    bool ok = (p[1] == 0.0 && p[4] == 0.0);
    std::printf("p after mode0 empty=(%g,%g) %s\n", p[1], p[4], ok?"PASS":"FAIL");
    return ok?0:1;
}
