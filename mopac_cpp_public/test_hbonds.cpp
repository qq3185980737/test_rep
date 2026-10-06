// test_hbonds.cpp — numeric check for MOZYME hydrogen-bond identification.
#include "hbonds.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"

#include <cstdio>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace MOZYME_C;

// local ijbo stub: atom pair (2,1) maps to packed block start 0 (i.e. first
// element of p). All other pairs are "not bonded" (-1).
int ijbo(int i, int j) { return (i == 2 && j == 1) ? 0 : -1; }

// shared storage for module-level MOZYME vectors used here (header declares
// them extern; single definition in this test binary)
namespace MOZYME_C { std::vector<int> iorbs; }
std::vector<int> ncf, nnce, nce, nncf, icocc, icvir;
std::vector<double> fmo;
std::vector<std::vector<int>> ifmo;

int main() {
    int ok = 1;
    numat = 2;
    norbs = 2;
    mpack = 1;
    iorbs.assign(3, 0); iorbs[1] = 1; iorbs[2] = 1;

    // density: empty so that sump < 1e-10 holds
    p.assign(2, 0.0);
    // Fock packed block: ll starts at ijbo(2,1)=0, then ++ll -> 1
    double fao_storage[2] = {0.0, 0.5};
    const double* fao = fao_storage;   // fao[1] = 0.5 -> sumf = 0.25

    // one occupied LMO containing atom 1 only (ncf=1 -> kk=0)
    int nocc = 1, nvir = 1;
    ncf.assign(2, 0); ncf[1] = 1;
    nncf.assign(2, 0); nncf[1] = 0;
    icocc.assign(3, 0); icocc[1] = 1;              // LMO 1 contains atom 1
    // one virtual LMO containing atom 2 only (nce=1 -> kk=0)
    nce.assign(2, 0); nce[1] = 1;
    nnce.assign(2, 0); nnce[1] = 0;
    icvir.assign(3, 0); icvir[1] = 2;              // virtual LMO 1 contains atom 2

    fmo.assign(4, 0.0);
    ifmo.assign(3, std::vector<int>(4, 0));

    int iused[8] = {0};
    int nij_loc = 0;
    hbonds(fao, nocc, nvir, iused, nij_loc, 0.01);

    if (nij_loc != 1) { std::printf("FAIL nij_loc=%d\n", nij_loc); ok = 0; }
    else {
        if (fmo[1] != 0.1) { std::printf("FAIL fmo[1]=%g\n", fmo[1]); ok = 0; }
        if (ifmo[2][1] != 1 || ifmo[1][1] != 1) { std::printf("FAIL ifmo: %d %d\n", ifmo[1][1], ifmo[2][1]); ok = 0; }
    }
    // cutoff too high -> nothing found
    nij_loc = 0;
    hbonds(fao, nocc, nvir, iused, nij_loc, 0.5);
    if (nij_loc != 0) { std::printf("FAIL cutoff: nij_loc=%d\n", nij_loc); ok = 0; }
    // populated density -> sump >= 1e-10 -> nothing found
    p.assign(2, 0.0); p[1] = 1e-3;
    nij_loc = 0;
    hbonds(fao, nocc, nvir, iused, nij_loc, 0.01);
    if (nij_loc != 0) { std::printf("FAIL density: nij_loc=%d\n", nij_loc); ok = 0; }

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
