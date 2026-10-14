// test_charst.cpp
#include <cstdio>
#include <vector>

#include "charst.h"

int main() {
    std::vector<std::vector<double>> v(1, std::vector<double>(1, 1.0));
    std::vector<int> nt = {0};
    std::vector<std::vector<double>> r(4, std::vector<double>(4, 0.0));
    bool first = true;

    double a = charst(v, nt, 1, 1, r, 1, first);
    bool r1 = (a == 1.0);
    std::printf("charst(E)=%.1f %s\n", a, r1 ? "PASS" : "FAIL");

    double b = charst(v, nt, -1, 1, r, 1, first);
    bool r2 = (b == 0.0);
    std::printf("charst(reset)=%.1f %s\n", b, r2 ? "PASS" : "FAIL");
    return (r1 && r2) ? 0 : 1;
}
