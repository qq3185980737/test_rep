// test_aabbcd.cpp — standalone unit test for aabbcd().
//
// Two hand-computed cases: case 2 also exercises the module side effect
// (ispqr write + is increment).
// Build with:
//   cl /std:c++17 /utf-8 /EHsc test_aabbcd.cpp aabbcd.cpp meci_C.cpp
// Exit code 0 == all PASS.

#include <cmath>
#include <cstdio>
#include <tuple>
#include <vector>

#include "aabbcd.h"
#include "meci_C.h"

namespace {

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

void check(const char* name, bool ok) {
    std::printf("%-36s %s\n", name, ok ? "PASS" : "FAIL");
    if (!ok) {
        ++failures;
    }
}

}  // namespace

int main() {
    // ---- Case 1: no side effect, even/odd phase computed manually. ----
    // iocca1 = {1,1,0,0}, iocca2 = {0,0,1,1}
    // ioccb1 = {1,0,1,0}, ioccb2 = {0,1,0,1}
    // i=0, j=1, k=0, l=1. Side effect: iocca1[0]==ioccb1[0] (1==1) -> skip.
    // No ordering swaps. xr = xy(0,1,0,1) = 5.
    // (i<=k && j<=l) -> ij=0; no additions.
    // sum over p=0..0: ioccb1[0]+iocca1[0] = 1+1 = 2; ij=2
    // sum over p=1..1: ioccb2[1]+iocca2[1] = 1+0 = 1; ij=3 (odd) -> -5
    {
        const std::vector<int> iocca1 = {1, 1, 0, 0};
        const std::vector<int> ioccb1 = {1, 0, 1, 0};
        const std::vector<int> iocca2 = {0, 0, 1, 1};
        const std::vector<int> ioccb2 = {0, 1, 0, 1};
        const auto xy = make_xy(4, {{0, 1, 0, 1, 5.0}});

        // Isolate module state for this case.
        meci_C::ispqr.assign(3, std::vector<int>(3, -99999));
        meci_C::iiloop = 1;
        meci_C::is = 1;
        meci_C::jloop = 42;

        const double got = aabbcd(iocca1.data(), ioccb1.data(),
                                  iocca2.data(), ioccb2.data(), 4, xy.data());
        check("case1 value == -5", std::fabs(got + 5.0) < 1e-12);
        check("case1 ispqr untouched (-99999)",
              meci_C::ispqr[1][1] == -99999 && meci_C::is == 1);
    }

    // ---- Case 2: side effect triggered, beta ordering swap exercised. ----
    // iocca1 = {1,0,1,1}, iocca2 = {0,1,1,1}
    // ioccb1 = {0,1,1,0}, ioccb2 = {1,0,1,0}
    // i=0, j=1, k=0, l=1. Side effect: i==k, j==l, iocca1[0]!=ioccb1[0]
    //   (1!=0) -> ispqr[1][1] = 42, is = 2.
    // alpha: no swap. beta: ioccb1[0]=0 < ioccb2[0]=1 -> swap k,l -> k=1,l=0.
    // xr = xy(0,1,1,0) = 6.
    // (i>k&&j>l)=false, (i<=k&&j<=l)=(0<=1 && 1<=0)=false -> ij=1
    // j>l (1>0): ij += iocca2[0]+ioccb2[1] = 0+0 = 0 -> ij=1
    // i>k? no. sum p=0..1: (0+1)+(1+0)=2 -> ij=3
    // j>l? yes -> swap j=0,l=1. sum p=0..1: (1+0)+(0+1)=2 -> ij=5 (odd) -> -6
    {
        const std::vector<int> iocca1 = {1, 0, 1, 1};
        const std::vector<int> ioccb1 = {0, 1, 1, 0};
        const std::vector<int> iocca2 = {0, 1, 1, 1};
        const std::vector<int> ioccb2 = {1, 0, 1, 0};
        const auto xy = make_xy(4, {{0, 1, 1, 0, 6.0}});

        meci_C::ispqr.assign(3, std::vector<int>(3, -99999));
        meci_C::iiloop = 1;
        meci_C::is = 1;
        meci_C::jloop = 42;

        const double got = aabbcd(iocca1.data(), ioccb1.data(),
                                  iocca2.data(), ioccb2.data(), 4, xy.data());
        check("case2 value == -6", std::fabs(got + 6.0) < 1e-12);
        check("case2 ispqr[1][1] == 42", meci_C::ispqr[1][1] == 42);
        check("case2 is incremented to 2", meci_C::is == 2);
    }

    if (failures == 0) {
        std::printf("\nAll checks PASSED.\n");
        return 0;
    }
    std::printf("\n%d check(s) FAILED.\n", failures);
    return 1;
}
