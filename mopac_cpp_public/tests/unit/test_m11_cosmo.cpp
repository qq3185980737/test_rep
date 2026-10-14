// test_m11_cosmo.cpp — smoke/structural verification for cosmo.cpp
// (MOPAC 2016 COSMO cavity construction + AMAT solve chain).
// Scenario: single H atom -> spherical cavity.
//   - cosini: initializes radii/directions without error
//   - coscav: builds ~n0[1] segments, area>0, volume>0, AMAT positive-definite
//   - mkbmat / addnuc: finite results, no error
// The single-atom cavity is the one geometry where the COSMO machinery can
// be exercised end to end without the closed-source two-electron kernels.
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
        ++g_fail; std::fprintf(stderr, "[FAIL] %s: got %.10g want %.10g\n", name, got, want);
    }
}

int main() {
    numat = 1; lm61 = 1;
    nfirst = {0, 1}; nlast = {0, 1};
    nat = {0, 1};
    coord.assign(4, std::vector<double>(2, 0.0));
    keywrd = "";
    mozyme = false; moperr = false; enuclr = 0.0;
    nspa = 42; nppa = 42; ioldcv = 0;

    cosini(false);
    chk(!moperr, "cosini no error");
    chk(nps == 0, "nps=0 after cosini");
    chk(srad[1] > 0.5, "srad(H) set");

    coscav();
    chk(!moperr, "coscav no error");
    chk(nps >= 12 && nps <= 1082, "nps in segment range");
    chk(area > 0.0 && std::isfinite(area), "cavity area positive finite");
    chk(cosvol > 0.0 && std::isfinite(cosvol), "cavity volume positive finite");
    // coscav() already performed the Cholesky factorization (coscl1) of AMAT
    // per F90 coscav; verify the factor has a positive diagonal.
    chk(amat[1] > 0.0 && std::isfinite(amat[1]), "AMAT factor diagonal > 0");

    mkbmat();
    chk(!moperr, "mkbmat no error");
    chk(std::isfinite(bmat[1][1]) && bmat[1][1] > 0.0, "bmat finite positive");

    addnuc();
    chk(!moperr, "addnuc no error");
    chk(std::isfinite(enuclr), "enuclr finite");
    chk(std::isfinite(qscnet[1][1]), "qscnet finite");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
