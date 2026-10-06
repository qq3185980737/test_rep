// test_digit.cpp
#include <cstdio>
#include <cmath>
#include "digit.h"
int main() {
    double a = digit("  123.456 rest", 3);
    double b = digit("-78.9xyz", 1);
    double c = digit("+42", 1);
    bool ok = std::fabs(a-123.456)<1e-9 && std::fabs(b+78.9)<1e-9 && std::fabs(c-42)<1e-9;
    std::printf("digit %g %g %g %s\n", a,b,c, ok?"PASS":"FAIL");
    return ok?0:1;
}
