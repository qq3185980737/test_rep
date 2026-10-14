// test_buildf.cpp
#include <cstdio>
#include <vector>

#include "buildf.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    numat = 1; mpack = 2; id = 0;
    f = {0, 0, 0};
    h = {0, 10.0, 20.0};
    std::vector<double> partf = {0, 5.0, 7.0};

    buildf(f, partf, 1);
    bool r1 = f[1] == 15.0 && f[2] == 27.0;
    std::printf("mode=1 f=(%.0f,%.0f) %s\n", f[1], f[2], r1 ? "PASS" : "FAIL");

    buildf(f, partf, -1);
    bool r2 = f[1] == -5.0 && f[2] == -13.0;
    std::printf("mode=-1 f=(%.0f,%.0f) %s\n", f[1], f[2], r2 ? "PASS" : "FAIL");

    buildf(f, partf, 0);
    bool r3 = f[1] == 10.0 && f[2] == 20.0;
    std::printf("mode=0 f=(%.0f,%.0f) %s\n", f[1], f[2], r3 ? "PASS" : "FAIL");

    return (r1 && r2 && r3) ? 0 : 1;
}
