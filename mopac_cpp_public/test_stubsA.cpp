// test_stubsA.cpp — smoke test for ionout/modgra/setupg/values/getsym.
#include "setupg.h"
#include "modgra.h"
#include "ionout.h"
#include "values.h"
#include "getsym.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>

void build_res_start_etc();
namespace molkst_C {
extern int norbs, mpack, natoms, numat, maxtxt, nvar, id, nclose, nopen, nelecs;
}
namespace common_arrays_C {
extern std::vector<int> nat, nfirst, nlast, labels, na, nb, nc;
extern std::vector<double> grad, q;
extern std::vector<std::vector<int>> loc;
extern std::vector<std::vector<double>> coord;
extern std::vector<std::string> txtatm;
}
namespace parameters_C {
extern double zs[108], zp[108];
}
namespace MOZYME_C {
extern std::vector<int> iorbs, res_start, at_res;
extern bool lijbo;
extern std::vector<std::vector<int>> nijbo;
}
extern std::vector<int> nce, ncf, ncocc, ncvir, nnce, nncf, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int nvirtual, noccupied;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;
extern void nuchar(char*, int, double*, int&);
extern int ijbo(int, int);

int main() {
    using common_arrays_C::nat;
    using common_arrays_C::nfirst;
    using common_arrays_C::nlast;
    using common_arrays_C::coord;
    using common_arrays_C::grad;
    using common_arrays_C::loc;
    using common_arrays_C::txtatm;
    using common_arrays_C::labels;
    using common_arrays_C::na;
    using common_arrays_C::nb;
    using common_arrays_C::nc;
    using molkst_C::numat;
    using molkst_C::nvar;
    using molkst_C::maxtxt;
    using molkst_C::natoms;
    using molkst_C::norbs;
    using molkst_C::mpack;
    using parameters_C::zs;
    using parameters_C::zp;

    bool ok = true;
    std::printf("A0\n"); std::fflush(stdout);
    // setupg: 1 atom H (ni=1, nqn=1), 1 s orbital only.
    natoms = 1; numat = 1;
    nat.resize(2); nfirst.resize(2); nlast.resize(2);
    nat[1] = 1; nfirst[1] = 1; nlast[1] = 1;
    coord.assign(3, std::vector<double>(2, 0.0));
    MOZYME_C::lijbo = true;
    MOZYME_C::nijbo.assign(2, std::vector<int>(2, 0));
    MOZYME_C::nijbo[1][1] = 0;  // ijbo = address-1: block for atom 1 starts at 1
    zs[1] = 1.188; zp[1] = 1.0;
    std::printf("A1\n"); std::fflush(stdout);
    setupg();
    std::printf("A2\n"); std::fflush(stdout);
    // Check master table values.
    if (!(std::fabs(overlaps_C::allz[1][1][1] - 23.10303149) < 1e-8)) {
        std::printf("FAIL allz(1,1)\n"); ok = false;
    }
    if (!(std::fabs(overlaps_C::allc[6][2][2] - 0.1051855189) < 1e-9)) {
        std::printf("FAIL allc(6,2,2)\n"); ok = false;
    }
    // H: ccc(1,1) = allc(1,1,1)=9.163596280e-3; zzz(1,1) = allz(1,1,1)*zs^2.
    if (!(std::fabs(overlaps_C::ccc[1][1] - 0.009163596280) < 1e-12)) {
        std::printf("FAIL ccc(1,1)\n"); ok = false;
    }
    if (!(std::fabs(overlaps_C::zzz[1][1] - 23.10303149 * 1.188 * 1.188) < 1e-6)) {
        std::printf("FAIL zzz(1,1)\n"); ok = false;
    }

    // modgra with a trivial 1-atom system (no residues).
    nvar = 0;
    txtatm.assign(2, std::string(40, ' '));
    txtatm[1] = "ATOM      1  H   MOL     1";
    MOZYME_C::at_res.assign(2, 0);
    MOZYME_C::res_start.assign(2, 0);
    modgra();
    std::printf("modgra smoke OK\n");

    // ionout with one ion.
    maxtxt = 0;
    labels.assign(2, 1);
    na.assign(2, 0); nb.assign(2, 0); nc.assign(2, 0);
    int ions[4] = {1, 0, 0, 1};
    ionout(ions, 1, 1);
    ionout(ions, 3, 1);
    std::printf("ionout smoke OK\n");

    // values with tiny MOZYME setup: n=1 occupied, diagonal f=1, cocc=1.
    norbs = 1; mpack = 1;
    noccupied = 1; nvirtual = 0;
    icocc_dim = 1; cocc_dim = 1;
    icocc.assign(2, 1); ncocc.assign(2, 0); nncf.assign(2, 0); ncf.assign(2, 1);
    cocc.assign(2, 1.0);
    MOZYME_C::iorbs.assign(2, 1);
    common_arrays_C::f.assign(2, 1.0);
    common_arrays_C::eigs.assign(2, 0.0);
    values("OCCUPIED");
    if (!(std::fabs(common_arrays_C::eigs[1] - 1.0) < 1e-9)) {
        std::printf("FAIL values eigs[1]=%g\n", common_arrays_C::eigs[1]); ok = false;
    }
    values("BADTYPE");
    std::printf("values smoke OK\n");

    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
