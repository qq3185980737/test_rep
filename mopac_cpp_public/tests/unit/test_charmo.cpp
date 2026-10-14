// test_charmo.cpp
#include <cstdio>
#include <vector>

#include "charmo.h"
#include "molkst_C.h"

int main() {
    molkst_C::norbs = 1; molkst_C::numat = 1;
    std::vector<std::vector<double>> v(2, std::vector<double>(2, 1.0));
    std::vector<int> nt = {0, 100};
    std::vector<std::vector<double>> r(4, std::vector<double>(4, 0.0));
    bool first = true;
    double val = charmo(v, nt, 1, 1, r, 1, first);
    bool ok = (val == 1.0);
    std::printf("charmo(E)=%.1f %s\n", val, ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
