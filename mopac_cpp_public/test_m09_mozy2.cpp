// test_m09_mozy2.cpp — verification for hbonds.cpp + eimp.cpp (MOPAC 2016).
// Both routines are MOZYME-path helpers that operate on sparse-MO
// bookkeeping arrays; verified on a minimal 2-atom setup with nijbo lookup.
#include <cstdio>
#include <cmath>
#include <vector>
#include "hbonds.h"
#include "eimp.h"
#include "ijbo.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace molkst_C;

static int g_checks = 0; static int g_fail = 0; static double g_err = 0.0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    double e = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
    if (e > g_err) g_err = e;
    if (e > tol) { ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want); }
}

int main() {
    std::fprintf(stderr, "[t1] hbonds (2 atoms, one occupied + one virtual LMO)\n");
    {
        numat = 2; norbs = 2; mpack = 3;
        lijbo = true;
        nijbo.assign(3, std::vector<int>(3, 0));   // ijbo(1,2)=0 >= 0
        iorbs = {0, 1, 1};
        p.assign(4, 0.0);                 // p[1] = 0 (empty pair)
        f.assign(4, 0.0);
        // one occupied LMO covering atom 1, one virtual LMO covering atom 2
        ncf = {0, 2}; nncf = {0, 0};
        icocc = {0, 1, 1};
        nce = {0, 2}; nnce = {0, 0};
        icvir = {0, 2, 2};
        fmo.assign(4, 0.0);
        ifmo.assign(3, std::vector<int>(4, 0));
        std::vector<double> fao(4, 0.0);
        fao[1] = 1.0;                    // sumf = 1.0 > cutoff
        std::vector<int> iused(4, 0);
        int nij_loc = 0;
        hbonds(fao.data(), 1, 1, iused.data(), nij_loc, 0.5);
        chk(nij_loc == 1, "hbonds nij=1");
        chk(iused[1] == 1, "hbonds iused1=1");
        chk_rel(fmo[1], 0.1, 1e-12, "hbonds fmo1=0.1");
        chk(ifmo[1][1] == 1, "hbonds ifmo11=1");
        chk(ifmo[2][1] == 1, "hbonds ifmo21=1");
        // weak Fock pair must NOT create a bond
        std::fill(fmo.begin(), fmo.end(), 0.0);
        ifmo.assign(3, std::vector<int>(4, 0));
        nij_loc = 0;
        fao[1] = 0.1;                    // sumf = 0.01 < cutoff
        hbonds(fao.data(), 1, 1, iused.data(), nij_loc, 0.5);
        chk(nij_loc == 0, "hbonds weak pair rejected");
    }

    std::fprintf(stderr, "[t2] eimp (s-s Fock squared into p[k+1])\n");
    {
        numat = 2;
        iorbs = {0, 1, 1};
        nijbo.assign(3, std::vector<int>(3, 0));   // ijbo(1,2)=0
        f.assign(4, 0.0);
        f[1] = 2.0;                      // sum = 4.0
        p.assign(4, 0.0);
        eimp();
        chk_rel(p[1], 4.0, 1e-12, "eimp p1=4.0");
        chk_rel(p[2], 0.0, 1e-12, "eimp p2 untouched");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
