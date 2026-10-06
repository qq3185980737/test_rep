// test_charvi.cpp
#include <cstdio>
#include <vector>

#include "charvi.h"

int main() {
    std::vector<std::vector<double>> v(4, std::vector<double>(2, 1.0));
    std::vector<std::vector<double>> r(4, std::vector<double>(4, 0.0));
    double a = charvi(v, 1, 1, r, 3);
    bool ok = (a == 1.0);
    std::printf("charvi(E)=%.1f %s\n", a, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
