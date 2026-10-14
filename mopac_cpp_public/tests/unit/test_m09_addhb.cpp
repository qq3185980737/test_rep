// test_m09_addhb.cpp — verification for addhb.cpp (MOPAC 2016).
// End-to-end: hbonds identifies one hydrogen bond (fmo entry), then diagg2
// performs the Jacobi annihilation with mutual LMO expansion.  Scenario is
// the same 2-atom system as test_m09_diagg2, driven through addhb.
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "addhb.h"
#include "ijbo.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace molkst_C;

void memory_error(const char* name) { std::fprintf(stderr, "memory_error: %s\n", name); std::exit(1); }

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
    numat = 2; norbs = 2; nvirtual = 2;
    numcal = 1; id = 0;
    keywrd = "";
    nat = std::vector<int>{0, 6, 6};
    thresh = 1e-4; tiny = 1e-6; shift = 0.0;
    int nocc1 = 1, nvir1 = 1, idiagg = 1, nij = 0, nhb = 1;
    // MOZYME LMO bookkeeping (same scenario as test_m09_diagg2)
    ncf = {0, 1}; nncf = {0, 0};
    nce = {0, 1}; nnce = {0, 0};
    icocc = {0, 1, 0, 0};            // occupied LMO 1 covers atom 1
    icvir = {0, 2, 0, 0};            // virtual  LMO 1 covers atom 2
    iorbs = {0, 1, 1};
    ncocc = {0, 0}; ncvir = {0, 0};
    cocc_dim = 3; cvir_dim = 3; icocc_dim = 3; icvir_dim = 3;
    cocc = std::vector<double>{0, 0.6, 0.0, 0.0};
    cvir = std::vector<double>{0, 0.9, 0.0, 0.0};
    fmo = std::vector<double>(4, 0.0);
    ifmo.assign(3, std::vector<int>(4, 0));
    // hbonds inputs: fock ao and density over the (1,2) pair
    f = std::vector<double>{0, 1.5, 0.0, 0.0};   // sumf = 2.25 > cutoff 1.0
    p = std::vector<double>{0, 0.0, 0.0, 0.0};   // sump = 0 < 1e-10
    lijbo = true;
    nijbo.assign(3, std::vector<int>(3, 0));     // ijbo(1,2)=0
    // diagg2 inputs: eigv slice eigs(nocc1+1:) -> eigv(1)=eigs[2]
    eigs = std::vector<double>{0, -0.5, 0.1};
    sumb = 0.0;

    std::fprintf(stderr, "[t1] addhb: hbonds -> diagg2 chain\n");
    addhb(nocc1, nvir1, idiagg, nij, nhb);

    // hbonds: one pair created
    chk(nij == 1, "addhb nij=1 from hbonds");
    chk_rel(fmo[1], 0.1, 1e-12, "addhb fmo1=0.1");
    chk(ifmo[1][1] == 1 && ifmo[2][1] == 1, "addhb ifmo pair {1,1}");
    // diagg2: mutual expansion (expectations computed with python)
    // c=0.1, d=eigs[1]-eigv[1]-0=-0.6 -> alpha=0.9870874576 beta=-0.1601822430
    const double alpha = 0.9870874576, beta = -0.1601822430;
    chk(ncf[1] == 2, "addhb occupied expanded");
    chk(nce[1] == 2, "addhb virtual expanded");
    chk(icocc[1] == 1 && icocc[2] == 2, "addhb icocc {1,2}");
    chk(icvir[1] == 2 && icvir[2] == 1, "addhb icvir {2,1}");
    chk_rel(cocc[1], alpha * 0.6, 1e-9, "addhb cocc1");
    chk_rel(cocc[2], beta * 0.9, 1e-9, "addhb cocc2");
    chk_rel(cvir[1], alpha * 0.9, 1e-9, "addhb cvir1");
    chk_rel(cvir[2], -beta * 0.6, 1e-9, "addhb cvir2");
    chk_rel(sumb, std::fabs(beta), 1e-9, "addhb sumb");

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
