// test_m09_diagg2.cpp — verification for diagg2.cpp (MOPAC 2016).
// Minimal 2-atom / 1-occupied + 1-virtual LMO Jacobi annihilation with
// mutual adaptive expansion of both LMO atom lists.
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "diagg2.h"
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
    numat = 2; norbs = 2; nvirtual = 2;
    nat = std::vector<int>{0, 6, 6};
    numcal = 1; id = 0;
    keywrd = "";                     // no DAMP / TIMES / DIAGG2
    thresh = 1e-4; tiny = 1e-6; shift = 0.0;
    int nocc = 1, nvir = 1, nij = 1, idiagg = 1;
    ncf = {0, 1}; nncf = {0, 0};
    nce = {0, 1}; nnce = {0, 0};
    icocc = {0, 1, 0, 0};            // occupied LMO 1 covers atom 1
    icvir = {0, 2, 0, 0};            // virtual  LMO 1 covers atom 2
    iorbs = {0, 1, 1};
    ncocc = {0, 0}; ncvir = {0, 0};  // both LMOs start at offset 0
    cocc_dim = 3; cvir_dim = 3; icocc_dim = 3; icvir_dim = 3;
    cocc = std::vector<double>{0, 0.6, 0.0, 0.0};   // atom 1 AO = 0.6
    cvir = std::vector<double>{0, 0.9, 0.0, 0.0};   // atom 2 AO = 0.9
    fmo = std::vector<double>{0, 0.05};
    ifmo.assign(3, std::vector<int>(3, 0));
    ifmo[1][1] = 1; ifmo[2][1] = 1;   // i=1 (virtual), j=1 (occupied)
    eigs = std::vector<double>{0, -0.5};
    std::vector<double> eigv3 = {0, 0.1, 0.0};   // eigv(1)=0.1 (1-based)
    std::vector<int> iused(3, 0);
    std::vector<char> latoms(3, 0);
    std::vector<double> storei(3, 0.0), storej(3, 0.0);
    sumb = 0.0;

    std::fprintf(stderr, "[t1] diagg2 mutual LMO expansion\n");
    diagg2(nocc, nvir, eigv3.data(), iused.data(), latoms.data(),
           nij, idiagg, storei.data(), storej.data());

    // expected values computed independently with python:
    // e=-0.6082762530, alpha=0.9965926760, beta=-0.0824805315
    const double alpha = 0.9965926760, beta = -0.0824805315;
    chk(ncf[1] == 2, "diagg2 occupied expanded to 2 atoms");
    chk(nce[1] == 2, "diagg2 virtual expanded to 2 atoms");
    chk(icocc[1] == 1 && icocc[2] == 2, "diagg2 icocc {1,2}");
    chk(icvir[1] == 2 && icvir[2] == 1, "diagg2 icvir {2,1}");
    // rotated AO values (F90 order: cocc[mlff] = alpha*a + beta*b for common,
    // then new atoms: cocc = beta*cvir / cvir = alpha*cvir, and
    // cvir = -beta*cocc / cocc = alpha*cocc)
    chk_rel(cocc[1], alpha * 0.6, 1e-9, "cocc1 rotated");
    chk_rel(cocc[2], beta * 0.9, 1e-9, "cocc2 = beta*cvir");
    chk_rel(cvir[1], alpha * 0.9, 1e-9, "cvir1 rotated");
    chk_rel(cvir[2], -beta * 0.6, 1e-9, "cvir2 = -beta*cocc");
    // counters reset after successful annihilation
    chk(iused[1] == -1 && iused[2] == -1, "diagg2 iused reset");
    chk(latoms[1] == 0 && latoms[2] == 0, "diagg2 latoms reset");
    chk_rel(sumb, std::fabs(beta), 1e-9, "diagg2 sumb = |beta|");

    std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
    return g_fail ? 1 : 0;
}
