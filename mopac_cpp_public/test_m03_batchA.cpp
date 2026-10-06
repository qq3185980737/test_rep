// test_m03_batchA.cpp — M03 batch A: CI two-electron integral family + bfn
// aababc / aabacd / aabbcd / babbbc / babbcd (MECI matrix elements) + bfn (B integrals)
#include <cmath>
#include <cstdio>
#include <vector>
#include "aababc.h"
#include "aabacd.h"
#include "aabbcd.h"
#include "babbbc.h"
#include "babbcd.h"
#include "bfn.h"
#include "meci_C.h"
#include "overlaps_C.h"

static int passed = 0, failed = 0;
static void check(const char* name, bool ok) {
    if (ok) { ++passed; std::printf("  [PASS] %s\n", name); }
    else { ++failed; std::printf("  [FAIL] %s\n", name); }
}
static bool close(double a, double b, double tol = 1e-10) {
    return std::fabs(a - b) < tol;
}
static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

int main() {
    std::printf("M03 batch A tests (CI integrals + bfn)\n");
    const int NM = 3;
    meci_C::occa.assign(4, 0.0);
    meci_C::occa[1] = 1.0; meci_C::occa[2] = 1.0; meci_C::occa[3] = 0.0;

    std::vector<double> xy(NM * NM * NM * NM, 0.0);
    for (int i = 1; i <= NM; ++i)
        for (int j = 1; j <= NM; ++j)
            for (int k = 1; k <= NM; ++k)
                for (int l = 1; l <= NM; ++l)
                    xy[xyidx(i, j, k, l, NM)] = i + 10.0 * j + 100.0 * k + 1000.0 * l;

    // 1) aababc: one alpha diff. alpha1={1,1,0} alpha2={1,0,1}
    // i=2, j=3 (immediate exit, no ij add); ij=b1[2]=0 even; v = sum formula
    {
        int a1[4] = {0, 1, 1, 0}, b1[4] = {0, 1, 0, 0};
        int a2[4] = {0, 1, 0, 1}, b2[4] = {0, 1, 0, 0};
        double v = aababc(a1, b1, a2, NM, xy.data());
        double s = 0;
        for (int k = 1; k <= NM; ++k) {
            double xa = xy[xyidx(2, 3, k, k, NM)], xb = xy[xyidx(2, k, 3, k, NM)];
            s += (xa - xb) * (a1[k] - meci_C::occa[k]) + xa * (b1[k] - meci_C::occa[k]);
        }
        check("aababc one-alpha element", close(v, s));
    }

    // 2) aabacd: NM=4. alpha1={1,1,0,0} alpha2={0,0,1,1}
    // i=3 (first a1<a2), j=4 (second); k=1 (first a1>a2), l=2 (second)
    // ij = b2[3]+b1[1] = 0 -> even
    // v = xy(i,k,j,l)-xy(i,l,k,j) = xy(3,1,4,2)-xy(3,2,1,4)  [l 2nd, k 3rd, j 4th index]
    {
        const int NM4 = 4;
        std::vector<double> xy4(NM4 * NM4 * NM4 * NM4, 0.0);
        for (int i = 1; i <= NM4; ++i)
            for (int j = 1; j <= NM4; ++j)
                for (int k = 1; k <= NM4; ++k)
                    for (int l = 1; l <= NM4; ++l)
                        xy4[xyidx(i, j, k, l, NM4)] = i + 10.0 * j + 100.0 * k + 1000.0 * l;
        int a1[5] = {0, 1, 1, 0, 0}, a2[5] = {0, 0, 0, 1, 1};
        int b1[5] = {0, 0, 0, 0, 0}, b2[5] = {0, 0, 0, 0, 0};
        double v = aabacd(a1, b1, a2, b2, NM4, xy4.data());
        double exp = xy4[xyidx(3, 1, 4, 2, NM4)] - xy4[xyidx(3, 2, 1, 4, NM4)];
        check("aabacd two-alpha element", close(v, exp));
    }

    // 3) aabbcd: alpha diffs i=2,j=3; beta diffs k=1,l=2
    // a1={1,1,0}->a2={1,0,1} ; b1={1,0,0}->b2={0,1,0}
    // xr = xy(2,3,1,2); ij even -> no flip
    {
        int a1[4] = {0, 1, 1, 0}, a2[4] = {0, 1, 0, 1};
        int b1[4] = {0, 1, 0, 0}, b2[4] = {0, 0, 1, 0};
        meci_C::ispqr.assign(4, std::vector<int>(4, 0));
        meci_C::iiloop = 1; meci_C::is = 1; meci_C::jloop = 7;
        double v = aabbcd(a1, b1, a2, b2, NM, xy.data());
        bool ok = (meci_C::ispqr[1][1] == 0);
        ok = ok && close(v, xy[xyidx(2, 3, 1, 2, NM)]);
        check("aabbcd one-alpha-one-beta element", ok);
    }

    // 4) aabbcd same-index path: i==k==1, j==l==3, a1[1]==b1[1] -> no ispqr; ij even -> no flip
    {
        int a1[4] = {0, 1, 1, 0}, a2[4] = {0, 0, 1, 1};
        int b1[4] = {0, 1, 0, 0}, b2[4] = {0, 0, 0, 1};
        meci_C::ispqr.assign(4, std::vector<int>(4, 0));
        meci_C::iiloop = 2; meci_C::is = 1; meci_C::jloop = 9;
        double v = aabbcd(a1, b1, a2, b2, NM, xy.data());
        bool ok = (meci_C::ispqr[2][1] == 0) && close(v, xy[xyidx(1, 3, 1, 3, NM)]);
        check("aabbcd same-index path", ok);
    }

    // 5) babbbc: beta1={1,0,0}->beta2={0,0,1}; i=1, j=3; ij=a1[2]+b1[2]=1 -> odd -> flip
    {
        int a1[4] = {0, 1, 1, 0}, b1[4] = {0, 1, 0, 0}, b2[4] = {0, 0, 0, 1};
        double v = babbbc(a1, b1, b2, NM, xy.data());
        double s = 0;
        for (int k = 1; k <= NM; ++k) {
            double xa = xy[xyidx(1, 3, k, k, NM)], xb = xy[xyidx(1, k, 3, k, NM)];
            s += (xa - xb) * (b1[k] - meci_C::occa[k]) + xa * (a1[k] - meci_C::occa[k]);
        }
        check("babbbc one-beta element", close(v, -s));
    }

    // 6) babbcd: NM=4. beta1={1,1,0,0}->beta2={0,0,1,1}
    // i=3, j=4; k=1, l=2; ij = a2[4]+a1[2] = 0 -> even -> one=+1
    // v = (xy(i,k,j,l)-xy(i,l,k,j))*one = (xy(3,1,4,2)-xy(3,2,1,4))*1
    {
        const int NM4 = 4;
        std::vector<double> xy4(NM4 * NM4 * NM4 * NM4, 0.0);
        for (int i = 1; i <= NM4; ++i)
            for (int j = 1; j <= NM4; ++j)
                for (int k = 1; k <= NM4; ++k)
                    for (int l = 1; l <= NM4; ++l)
                        xy4[xyidx(i, j, k, l, NM4)] = i + 10.0 * j + 100.0 * k + 1000.0 * l;
        int a1[5] = {0, 0, 0, 0, 0}, b1[5] = {0, 1, 1, 0, 0};
        int a2[5] = {0, 0, 0, 0, 0}, b2[5] = {0, 0, 0, 1, 1};
        double v = babbcd(a1, b1, a2, b2, NM4, xy4.data());
        double exp = (xy4[xyidx(3, 1, 4, 2, NM4)] - xy4[xyidx(3, 2, 1, 4, NM4)]) * (1.0);
        check("babbcd two-beta element", close(v, exp));
    }

    // 7) bfn B integrals (1-based bf[1..13], bf[0] padding)
    {
        overlaps_C::fact[0] = 1.0;
        for (int i = 1; i <= 17; ++i) overlaps_C::fact[i] = overlaps_C::fact[i - 1] * i;
        double bf[14] = {0};
        bfn(0.0001, bf);  // series (last=6)
        bool ok = close(bf[1], 2.000000003, 1e-7) && close(bf[2], -6.666666673e-05, 1e-9);
        bfn(1.5, bf);     // series (last=12)
        ok = ok && close(bf[1], 2.839039273, 1e-7) && close(bf[2], -1.2438533, 1e-7);
        bfn(3.5, bf);     // recurrence
        ok = ok && close(bf[1], 9.452929879, 1e-7) && close(bf[2], -6.769348418, 1e-7);
        bfn(1e-7, bf);    // small-x (label 90)
        ok = ok && close(bf[1], 2.0, 1e-12) && close(bf[3], 0.6666666666666666, 1e-12);
        check("bfn three branches", ok);
    }

    std::printf("M03 batch A: %d passed, %d failed\n", passed, failed);
    return failed ? 1 : 0;
}
