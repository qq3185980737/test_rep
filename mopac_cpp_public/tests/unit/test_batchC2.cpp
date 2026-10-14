// test_batchC2.cpp — tests: setup_mopac_arrays (mode 1 essentials, mode 2 main
// arrays, n=0 teardown). meci/fock2/mopend are stubbed; setup_mopac_arrays real.
#include <cstdio>
#include <string>
#include <vector>
#include "setup_mopac_arrays.h"
#include "common_arrays_C.h"
#include "maps_C.h"
#include "iter_C.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "derivs_C.h"
#include "ef_C.h"
#include "meci_C.h"
#include "drc_C.h"
#include "cosmo_C.h"
#include "esp_C.h"
#include "to_screen_C.h"

namespace molkst_C {
extern int numat, norbs, mpack, n2elec, nvar, l123, natoms, num_bits;
extern bool uhf, mozyme;
}
namespace chanel_C { extern int iw; }

double meci() { return 0.0; }

int fock2_calls = 0;
void fock2(std::vector<double>&, const std::vector<double>&, std::vector<double>&,
           const std::vector<double>&, const std::vector<double>&,
           const std::vector<double>&, int, const std::vector<int>&,
           const std::vector<int>&, int) { ++fock2_calls; }

using namespace common_arrays_C;
using namespace molkst_C;
using namespace iter_C;
using namespace symmetry_C;
using namespace maps_C;
using namespace derivs_C;
using namespace ef_C;
using namespace meci_C;
using namespace drc_C;
using namespace cosmo_C;
using namespace esp_C;
using namespace to_screen_C;

static void reset_globals() {
    numat = 0; norbs = 0; mpack = 0; n2elec = 0; nvar = 0; l123 = 1;
    natoms = 0; uhf = false; mozyme = false; num_bits = 64;
    // clear everything first, as setup(0,0) would
    setup_mopac_arrays(0, 0);
    fock2_calls = 0;
}

int main() {
    bool ok = true;

    // ---- T1: mode=1 essential arrays sized and initialised.
    {
        reset_globals();
        setup_mopac_arrays(3, 1);
        if (geo.size() != 3 || geo[0].size() != 4) { std::fprintf(stderr, "FAIL T1 geo=%zux%zu\n", geo.size(), geo[0].size()); ok = false; }
        if (coord.size() != 3 || coord[0].size() != 4) { std::fprintf(stderr, "FAIL T1 coord\n"); ok = false; }
        if (nfirst.size() != 4 || nlast.size() != 4 || nat.size() != 4) { std::fprintf(stderr, "FAIL T1 nfirst/nlast/nat sizes\n"); ok = false; }
        if (uspd.size() != 9 * 3 + 1 || pdiag.size() != 9 * 3 + 1) { std::fprintf(stderr, "FAIL T1 uspd/pdiag\n"); ok = false; }
        if (xparam.size() != 3 * 3 + 1 || xparef.size() != 3 * 3 + 1) { std::fprintf(stderr, "FAIL T1 xparam/xparef\n"); ok = false; }
        if (loc.size() != 2 || loc[0].size() != 3 * 3 + 1) { std::fprintf(stderr, "FAIL T1 loc\n"); ok = false; }
        if (ibonds.size() != 16 || ibonds[0].size() != 4) { std::fprintf(stderr, "FAIL T1 ibonds\n"); ok = false; }
        if (jelem.size() != 21 || jelem[0].size() != 4) { std::fprintf(stderr, "FAIL T1 jelem\n"); ok = false; }
        if (nfirst[1] != -9999 || nlast[3] != -9999) { std::fprintf(stderr, "FAIL T1 nfirst/nlast sentinel\n"); ok = false; }
        if (na_store[2] != 0 || nbonds[1] != 0) { std::fprintf(stderr, "FAIL T1 zeros\n"); ok = false; }
    }
    // ---- T2: mode=2 conventional (non-MOZYME) main arrays sized + zeroed.
    {
        reset_globals();
        numat = 2; norbs = 5; mpack = 15; n2elec = 100; nvar = 6; l123 = 1; natoms = 2;
        setup_mopac_arrays(2, 2);
        if (h.size() != 16 || p.size() != 16 || pa.size() != 16 || pb.size() != 16) { std::fprintf(stderr, "FAIL T2 h/p/pa/pb\n"); ok = false; }
        if (pold.size() != 6 * 15 + 1 || pold2.size() != 6 * 15 + 1) { std::fprintf(stderr, "FAIL T2 pold\n"); ok = false; }
        if (f.size() != 16) { std::fprintf(stderr, "FAIL T2 f\n"); ok = false; }
        if (c.size() != 6 || c[0].size() != 6) { std::fprintf(stderr, "FAIL T2 c\n"); ok = false; }
        if (eigs.size() != 7) { std::fprintf(stderr, "FAIL T2 eigs=%zu\n", eigs.size()); ok = false; }
        if (q.size() != 3) { std::fprintf(stderr, "FAIL T2 q=%zu\n", q.size()); ok = false; }
        if (eigb.size() != 6) { std::fprintf(stderr, "FAIL T2 eigb=%zu\n", eigb.size()); ok = false; }
        if (pold3.size() != std::max(mpack, 400) + 1) { std::fprintf(stderr, "FAIL T2 pold3=%zu\n", pold3.size()); ok = false; }
        if (w.size() != 100 + 2025 + 1) { std::fprintf(stderr, "FAIL T2 w=%zu\n", w.size()); ok = false; }
        if (dxyz.size() != 3 * 2 * 1 + 1) { std::fprintf(stderr, "FAIL T2 dxyz=%zu\n", dxyz.size()); ok = false; }
        if (grad.size() != 7) { std::fprintf(stderr, "FAIL T2 grad=%zu\n", grad.size()); ok = false; }
        if (errfn.size() != 3 * 2 * 1 + 1) { std::fprintf(stderr, "FAIL T2 errfn=%zu\n", errfn.size()); ok = false; }
        if (pold[3] != 0.0 || eigb[2] != 0.0) { std::fprintf(stderr, "FAIL T2 zero init\n"); ok = false; }
        if (!wk.empty()) { std::fprintf(stderr, "FAIL T2 wk should be empty (l123=1)\n"); ok = false; }
        if (fb.size() != 0 || cb.size() != 0) { std::fprintf(stderr, "FAIL T2 fb/cb should be empty (rhf)\n"); ok = false; }
    }
    // ---- T3: mode=2 with l123=2 -> wk allocated; uhf -> fb/cb/pbold allocated.
    {
        reset_globals();
        numat = 2; norbs = 5; mpack = 15; n2elec = 100; nvar = 6; l123 = 2; natoms = 2; uhf = true;
        setup_mopac_arrays(2, 2);
        if (wk.size() != 100 + 2025 + 1) { std::fprintf(stderr, "FAIL T3 wk=%zu\n", wk.size()); ok = false; }
        if (fb.size() != 16 || cb.size() != 6) { std::fprintf(stderr, "FAIL T3 fb/cb\n"); ok = false; }
        if (pbold.size() != 6 * 15 + 1 || pbold2.size() != 6 * 15 + 1 || pbold3.size() != std::max(mpack, 400) + 1) { std::fprintf(stderr, "FAIL T3 pbold*\n"); ok = false; }
        if (pbold[5] != 0.0) { std::fprintf(stderr, "FAIL T3 pbold zero\n"); ok = false; }
    }
    // ---- T4: MOZYME job: mode=2 allocates grad only.
    {
        reset_globals();
        numat = 2; norbs = 5; mpack = 15; n2elec = 100; nvar = 6; natoms = 2; mozyme = true;
        h.assign(100, 1.0);  // mark
        setup_mopac_arrays(2, 2);
        if (grad.size() != 7) { std::fprintf(stderr, "FAIL T4 grad=%zu\n", grad.size()); ok = false; }
        if (h.size() != 100 || h[5] != 1.0) { std::fprintf(stderr, "FAIL T4 h should be untouched (mozyme)\n"); ok = false; }
    }
    // ---- T5: teardown n=0 clears everything, resets globals, calls fock2.
    {
        reset_globals();
        numat = 2; norbs = 5; mpack = 15; n2elec = 100; nvar = 6; natoms = 2;
        setup_mopac_arrays(2, 2);
        setup_mopac_arrays(0, 0);
        if (mpack != 0 || numat != 0 || n2elec != 0) { std::fprintf(stderr, "FAIL T5 globals\n"); ok = false; }
        if (!h.empty() || !w.empty() || !p.empty() || !grad.empty() || !f.empty()) { std::fprintf(stderr, "FAIL T5 arrays not cleared\n"); ok = false; }
        if (!geo.empty() || !nat.empty() || !uspd.empty()) { std::fprintf(stderr, "FAIL T5 essentials not cleared\n"); ok = false; }
        if (fock2_calls != 1) { std::fprintf(stderr, "FAIL T5 fock2 calls=%d\n", fock2_calls); ok = false; }
    }
    std::fprintf(stderr, ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
