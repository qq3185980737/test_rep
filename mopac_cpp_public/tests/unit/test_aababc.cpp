// test_aababc.cpp — standalone unit test for aababc().
//
// The three cases below were hand-computed from the Fortran algorithm, so
// they do not depend on a Fortran compiler. Build with:
//   g++ -std=c++17 -Wall -Wextra -o test_aababc test_aababc.cpp aababc.cpp meci_C.cpp
// and run ./test_aababc (exit code 0 == all PASS).

#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

#include "aababc.h"
#include "meci_C.h"

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
    // ---- Case 1: nmos = 2, first difference at i=0, j=1, odd phase. ----
    // iocca1 = {1,0}, iocca2 = {0,1}, ioccb1 = {1,1}, occa = {2,1}
    // xy(0,1,0,0)=3, xy(0,0,1,0)=5; rest zero.
    // k=0: (3-5)*(1-2) + 3*(1-2) = 2-3 = -1
    // k=1: (0-0)*(0-1) + 0*(1-1) = 0
    // sum = -1, ij = ioccb1(0) = 1 (odd) -> sign flip -> +1
    {
        const std::vector<int> iocca1 = {1, 0};
        const std::vector<int> ioccb1 = {1, 1};
        const std::vector<int> iocca2 = {0, 1};
        meci_C::occa = {2.0, 1.0};
        const auto xy = make_xy(2, {{0, 1, 0, 0, 3.0}, {0, 0, 1, 0, 5.0}});
        check("case1 (odd phase)", aababc(iocca1.data(), ioccb1.data(),
                                          iocca2.data(), 2, xy.data()),
              1.0);
    }

    // ---- Case 2: nmos = 2, first difference at i=0, j=1, even phase. ----
    // iocca1 = {0,1}, iocca2 = {1,0}, ioccb1 = {0,1}, occa = {2,1}
    // k=0: (3-5)*(0-2) + 3*(0-2) = 4-6 = -2
    // k=1: 0
    // sum = -2, ij = ioccb1(0) = 0 (even) -> no flip -> -2
    {
        const std::vector<int> iocca1 = {0, 1};
        const std::vector<int> ioccb1 = {0, 1};
        const std::vector<int> iocca2 = {1, 0};
        meci_C::occa = {2.0, 1.0};
        const auto xy = make_xy(2, {{0, 1, 0, 0, 3.0}, {0, 0, 1, 0, 5.0}});
        check("case2 (even phase)", aababc(iocca1.data(), ioccb1.data(),
                                           iocca2.data(), 2, xy.data()),
              -2.0);
    }

    // ---- Case 3: nmos = 3, first difference at i=1 (i=0 agrees), j=2. ----
    // iocca1 = {1,1,0}, iocca2 = {1,0,1}, ioccb1 = {1,1,1}, occa = {2,2,1}
    // xy(1,2,0,0)=1, xy(1,0,2,0)=3, xy(1,2,1,1)=2, xy(1,1,2,1)=4
    // k=0: (1-3)*(1-2) + 1*(1-2) = 2-1 = 1
    // k=1: (2-4)*(1-2) + 2*(1-2) = 2-2 = 0
    // k=2: (0-0)*(0-1) + 0*(1-1) = 0
    // sum = 1, ij = ioccb1(1) = 1 (odd) -> sign flip -> -1
    {
        const std::vector<int> iocca1 = {1, 1, 0};
        const std::vector<int> ioccb1 = {1, 1, 1};
        const std::vector<int> iocca2 = {1, 0, 1};
        meci_C::occa = {2.0, 2.0, 1.0};
        const auto xy = make_xy(3, {{1, 2, 0, 0, 1.0},
                                    {1, 0, 2, 0, 3.0},
                                    {1, 2, 1, 1, 2.0},
                                    {1, 1, 2, 1, 4.0}});
        check("case3 (skipped first MO)",
              aababc(iocca1.data(), ioccb1.data(), iocca2.data(), 3, xy.data()),
              -1.0);
    }

    if (failures == 0) {
        std::printf("\nAll %d checks PASSED.\n", 3);
        return 0;
    }
    std::printf("\n%d check(s) FAILED.\n", failures);
    return 1;
}
