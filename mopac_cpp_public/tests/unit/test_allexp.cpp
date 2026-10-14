// test_allexp.cpp — for each x, print 4 exp variants; python compares vs FTN.
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <mkl.h>
extern "C" double __svml_exp1(double);
static void hex(const char* tag, double v) {
    std::uint64_t u; std::memcpy(&u, &v, 8);
    std::printf("%s %016llx\n", tag, (unsigned long long)u);
}
int main(int argc, char** argv) {
    // stdin: one hex x per line; stdout: 4 variants
    char buf[64];
    while (std::fgets(buf, sizeof buf, stdin)) {
        if (buf[0] == '\n' || buf[0] == 0) continue;
        std::uint64_t u = (std::uint64_t)std::strtoull(buf, nullptr, 16);
        double x; std::memcpy(&x, &u, 8);
        double na = -x;
        std::printf("X %016llx\n", (unsigned long long)u);
        hex("E_STD_N", std::exp(na));
        hex("E_SVM_N", __svml_exp1(na));
        hex("E_RSTD", 1.0 / std::exp(x));
        hex("E_RSVM", 1.0 / __svml_exp1(x));
    }
    return 0;
}
