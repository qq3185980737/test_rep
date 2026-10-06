// test_aabacd.cpp — standalone unit test for aabacd().
//
// The three cases below were hand-computed from the Fortran algorithm.
// Build with:
//   cl /std:c++17 /EHsc test_aabacd.cpp aabacd.cpp
// (or g++ -std=c++17 -Wall -o test_aabacd test_aabacd.cpp aabacd.cpp)
// Exit code 0 == all PASS.

#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

#include "aabacd.h"

namespace {

// Build a flat xy array of nmos^4 doubles in Fortran (column-major) order and
// set the listed 0-based (i,j,k,l) positions to the given values.
std::vector<double> make_xy(
    int nmos,
    std::initializer_list<std::tuple<int, int, int, int, double>> specs) {
    std::vector<double> xy(nmos * nmos * nmos * nmos, 0.0);
    for (const auto& s : specs) {
        const int i = std::get<0>(s);
        const int j = std::get<1>(s);
        const int k = std::get<2>(s);
        const int l = std::get<3>(s);
        const double v = std::get<4>(s);
        xy[i + nmos * (j + nmos * (k + nmos * l))] = v;
    }
    return xy;
}

int failures = 0;

void check(const char* name, double got, double expected) {
    const bool ok = std::fabs(got - expected) < 1e-12;
    std::printf("%-32s got=%+.12f  expected=%+.12f  %s\n", name, got, expected,
                ok ? "PASS" : "FAIL");
    if (!ok) {
        ++failures;
    }
}

}  // namespace

int main() {
    // ---- Case 1: nmos = 4, even phase. ----
    // iocca1 = {1,0,1,0}, iocca2 = {0,1,0,1}
    // ioccb1 = {1,1,0,1}, ioccb2 = {0,1,1,1}
    // i=1 (first iocca1<iocca2); j=3; ij after j-loop = iocca2[2]+ioccb2[2] = 1
    // k=0 (first iocca1>iocca2); l=2; ij after l-loop = +iocca1[1]+ioccb1[1] = 1
    // ij += ioccb2[1] + ioccb1[0] = 1+1 = 2  -> ij = 4 (even)
    // sum = xy(1,0,3,2) - xy(1,2,0,3) = 7 - 3 = 4
    {
        const std::vector<int> iocca1 = {1, 0, 1, 0};
        const std::vector<int> ioccb1 = {1, 1, 0, 1};
        const std::vector<int> iocca2 = {0, 1, 0, 1};
        const std::vector<int> ioccb2 = {0, 1, 1, 1};
        const auto xy = make_xy(4, {{1, 0, 3, 2, 7.0}, {1, 2, 0, 3, 3.0}});
        check("case1 (even phase)",
              aabacd(iocca1.data(), ioccb1.data(), iocca2.data(),
                     ioccb2.data(), 4, xy.data()),
              4.0);
    }

    // ---- Case 2: nmos = 4, odd phase. ----
    // Same as case 1 but ioccb1[0] = 0, so ij = 3 (odd) -> sign flip.
    {
        const std::vector<int> iocca1 = {1, 0, 1, 0};
        const std::vector<int> ioccb1 = {0, 1, 0, 1};
        const std::vector<int> iocca2 = {0, 1, 0, 1};
        const std::vector<int> ioccb2 = {0, 1, 1, 1};
        const auto xy = make_xy(4, {{1, 0, 3, 2, 7.0}, {1, 2, 0, 3, 3.0}});
        check("case2 (odd phase)",
              aabacd(iocca1.data(), ioccb1.data(), iocca2.data(),
                     ioccb2.data(), 4, xy.data()),
              -4.0);
    }

    // ---- Case 3: nmos = 5, immediate differences, no middle accumulation. ----
    // iocca1 = {1,1,0,0,0}, iocca2 = {0,0,1,1,0}
    // ioccb1 = {1,0,1,1,1}, ioccb2 = {1,1,0,1,0}
    // i=2, j=3; k=0, l=1 (j/l loops exit immediately, ij=0)
    // ij += ioccb2[2] + ioccb1[0] = 0+1 = 1 (odd)
    // sum = xy(2,0,3,1) - xy(2,1,0,3) = 5 - 2 = 3  -> -3
    {
        const std::vector<int> iocca1 = {1, 1, 0, 0, 0};
        const std::vector<int> ioccb1 = {1, 0, 1, 1, 1};
        const std::vector<int> iocca2 = {0, 0, 1, 1, 0};
        const std::vector<int> ioccb2 = {1, 1, 0, 1, 0};
        const auto xy = make_xy(5, {{2, 0, 3, 1, 5.0}, {2, 1, 0, 3, 2.0}});
        check("case3 (nmos=5, no middle accum.)",
              aabacd(iocca1.data(), ioccb1.data(), iocca2.data(),
                     ioccb2.data(), 5, xy.data()),
              -3.0);
    }

    if (failures == 0) {
        std::printf("\nAll %d checks PASSED.\n", 3);
        return 0;
    }
    std::printf("\n%d check(s) FAILED.\n", failures);
    return 1;
}
