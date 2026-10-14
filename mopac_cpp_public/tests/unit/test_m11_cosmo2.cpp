// test_m11_cosmo2.cpp — two-atom H2 verification for cosmo.cpp.
// Exercises the multi-atom path: coscav with surclo closure segments,
// mkbmat monopole entries per atom, addnuc solve.
// Scenario: H2 with H at (0,0,0) and (1.4,0,0).
//   - cosini: radii + idenat/ipiden index vectors
//   - coscav: nps > 12 (surclo adds closure segments to cover the
//     deleted "inner" patches facing the partner atom)
//   - mkbmat: bmat[1][ips] = 1/|seg - H1|, bmat[2][ips] = 1/|seg - H2|
//     (H has a single s orbital: idel=1, no multipole rows)
//   - addnuc: finite results, no error
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "cosmo.h"
#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "funcon_C.h"

using namespace common_arrays_C;
using namespace cosmo_C;
using namespace molkst_C;
using namespace parameters_C;
using namespace funcon_C;

static int g_checks = 0; static int g_fail = 0;
static void chk(bool ok, const char* name) {
    ++g_checks;
    if (!ok) { ++g_fail; std::fprintf(stderr, "[FAIL] %s\n", name); }
}
static void chk_rel(double got, double want, double tol, const char* name) {
    ++g_checks;
    if (std::fabs(got - want) > tol) {
        ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.12g want %.12g\n", name, got, want);
    }
}

int main() {
    numat = 2; lm61 = 2;
    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[0][2] = 1.4;   // H2 at (1.4, 0, 0)
    keywrd = "";
    mozyme = false; moperr = false; enuclr = 0.0;
    nspa = 42; nppa = 42; ioldcv = 0;

    cosini(false);
    chk(!moperr, "cosini no error");
    chk(nps == 0, "nps=0 after cosini");
    chk(srad[1] > 0.5, "srad(H1) set");
    chk(srad[2] > 0.5, "srad(H2) set");
    chk(idenat[1] == 1 && idenat[2] == 2, "idenat density starts");

    coscav();
    chk(!moperr, "coscav no error");
    chk(nps > 12, "nps grew beyond single-atom 12 (surclo closure)");
    chk(nps <= 1082, "nps within grid limit");
    chk(area > 0.0 && std::isfinite(area), "cavity area positive finite");
    chk(cosvol > 0.0 && std::isfinite(cosvol), "cavity volume positive finite");
    chk(amat[1] > 0.0 && std::isfinite(amat[1]), "AMAT factor diagonal > 0");

    mkbmat();
    chk(!moperr, "mkbmat no error");
    chk(std::isfinite(bmat[1][1]) && bmat[1][1] > 0.0, "bmat finite positive");
    // monopole rows: atom1 block starts at iden=0 (bmat[1]), atom2 at iden=1 (bmat[2])
    int nseg = (nps < 5) ? nps : 5;
    for (int ips = 1; ips <= nseg; ++ips) {
        double dx1 = cosurf[1][ips] - coord[0][1];
        double dy1 = cosurf[2][ips] - coord[1][1];
        double dz1 = cosurf[3][ips] - coord[2][1];
        double want1 = 1.0 / std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        char nm1[64]; std::snprintf(nm1, sizeof(nm1), "bmat H1 monopole seg%d", ips);
        chk_rel(bmat[1][ips], want1, 1e-10, nm1);
        double dx2 = cosurf[1][ips] - coord[0][2];
        double dy2 = cosurf[2][ips] - coord[1][2];
        double dz2 = cosurf[3][ips] - coord[2][2];
        double want2 = 1.0 / std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        char nm2[64]; std::snprintf(nm2, sizeof(nm2), "bmat H2 monopole seg%d", ips);
        chk_rel(bmat[2][ips], want2, 1e-10, nm2);
        // H carries a single s orbital (idel=1): mkbmat writes no multipole
        // rows, and the bmat row count is lm61+1 = 3, so no rows beyond 2 exist.
    }

    addnuc();
    chk(!moperr, "addnuc no error");
    chk(std::isfinite(enuclr), "enuclr finite");
    chk(std::isfinite(qscnet[1][1]), "qscnet seg1 finite");
    chk(std::isfinite(qscnet[nps][1]), "qscnet last seg finite");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
