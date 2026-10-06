// test_densit.cpp
#include <cstdio>
#include <vector>
#include "densit.h"
int main() {
    std::vector<std::vector<double>> c(3, std::vector<double>(3, 0.0));
    c[1][1] = 1.0; c[2][2] = 1.0;
    std::vector<double> p(4, 0.0);
    densit(c, 3, 2, 1, 2.0, 1, 0.0, p, 0);    bool ok = (p[1] == 2.0 && p[2] == 0.0 && p[3] == 0.0);
    std::printf("p=(%g,%g,%g) %s\n", p[1], p[2], p[3], ok?"PASS":"FAIL");
    return ok?0:1;
}
