// test_perm.cpp — numerical verification of perm() microstate generator.
// perm(iperm, nels, nmos, nperms, limci) enumerates the C(nmos,nels)
// combinations as columns of 0/1 occupancy strings.
#include <cstdio>
#include <cmath>
#include <vector>
#include "perm.h"

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}

int main() {
    // Case 1: nels=2, nmos=4 -> 6 permutations.
    {
        int nmos = 4, nels = 2;
        std::vector<int> ip(nmos * 128, 0);
        int nperms = 128;   // storage cap (F90: maxci*4)
        perm(ip.data(), nels, nmos, nperms, 0);
        chk(nperms == 6, "C(4,2)=6 perms");
        int count[4] = {0,0,0,0};
        bool dup = false;
        for (int p = 1; p <= nperms; ++p) {
            int ones = 0;
            for (int i = 1; i <= nmos; ++i) if (ip[(p-1)*nmos + (i-1)]) ++ones;
            if (ones != nels) chk(false, "each perm has nels ones");
            for (int i = 1; i <= nmos; ++i) count[i-1] += ip[(p-1)*nmos + (i-1)];
            for (int q = p + 1; q <= nperms; ++q) {
                bool same = true;
                for (int i = 1; i <= nmos; ++i)
                    if (ip[(p-1)*nmos+(i-1)] != ip[(q-1)*nmos+(i-1)]) same = false;
                if (same) dup = true;
            }
        }
        chk(!dup, "no duplicate columns");
        for (int i = 0; i < nmos; ++i) chk(count[i] == 3, "each MO occupied in 3 perms");
        chk(true, "C(4,2) structure ok");
    }
    // Case 2: nels=3, nmos=6 -> 20.
    {
        int nmos = 6, nels = 3;
        std::vector<int> ip(nmos * 512, 0);
        int nperms = 512;
        perm(ip.data(), nels, nmos, nperms, 0);
        chk(nperms == 20, "C(6,3)=20 perms");
    }
    // Case 3: limci=2 on nels=2, nmos=4.
    // Ground state 1100; single excitations 1010,1001,0110,0101 have dist 2;
    // 0011 has dist 4.  limci=2 keeps 5.
    {
        int nmos = 4, nels = 2;
        std::vector<int> ip(nmos * 128, 0);
        int nperms = 128;
        perm(ip.data(), nels, nmos, nperms, 2);
        chk(nperms == 5, "limci=2 keeps 5 microstates");
        for (int p = 1; p <= nperms; ++p) {
            int d = 0;
            for (int j = 1; j <= nmos; ++j) d += std::abs(ip[(p-1)*nmos+(j-1)] - ip[(0)*nmos+(j-1)]);
            chk(d <= 2, "all kept states within limci of ground state");
        }
    }
    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
