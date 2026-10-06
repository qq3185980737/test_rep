// test_linear_cosmo_batchA.cpp — verify linear_cosmo batch A:
// ini_linear_cosmo (allocation + iatom_pos/ijbo_diag), bpnew_vec (s-only
// density path), bz_vec (core-charge interaction), get_bvec (multipoles),
// mult_triangle_vec (packed triangle multiply), some_norm.
#include <cmath>
#include <cstdio>
#include <vector>
#include "linear_cosmo.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "MOZYME_C.h"

using namespace molkst_C;
using namespace cosmo_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;
using namespace MOZYME_C;
using namespace linear_cosmo;

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
    // --- setup: H2, single s orbital per atom ---
    numat = 2; norbs = 2; lenabc = 200; nps = 2;
    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[0][1] = 0.0; coord[1][1] = 0.0; coord[2][1] = 0.0;
    coord[0][2] = 1.4; coord[1][2] = 0.0; coord[2][2] = 0.0;
    tore[1] = 1.0; dd[1] = 1.0; qq[1] = 1.0;
    // ijbo via MOZYME lookup: zero offsets
    lijbo = true;
    nijbo.assign(3, std::vector<int>(3, 0));

    // --- ini_linear_cosmo ---
    ini_linear_cosmo();
    chk(!rsc.empty() && rsc.size() == (size_t)(4 * 70 * numat), "rsc allocated (4*70*numat)");
    chk(linear_cosmo::tm.size() == (size_t)(9 * numat), "tm allocated (3*3*numat)");
    chk(iatom_pos[1] == 0, "iatom_pos(1)=0");
    chk(iatom_pos[2] == 1, "iatom_pos(2)=1 (single orbital on atom 1)");
    chk(nipsrs.size() == (size_t)(lenabc + 1), "nipsrs sized lenabc");
    chk(solv_energy == 0.0, "solv_energy reset");

    // --- bz_vec: cores vs induced charges ---
    cosurf.assign(5, std::vector<double>(3, 0.0));
    cosurf[1][1] = 0.5; cosurf[2][1] = 0.0; cosurf[3][1] = 0.0;
    cosurf[1][2] = 2.0; cosurf[2][2] = 0.0; cosurf[3][2] = 0.0;
    std::vector<double> v(3, 0.0);
    bz_vec(v);
    chk_rel(v[1], 1.0 / 0.5 + 1.0 / 0.9, 1e-10, "bz_vec v1 = sum tore/r");
    chk_rel(v[2], 1.0 / 2.0 + 1.0 / 0.6, 1e-10, "bz_vec v2 = sum tore/r");

    // --- bpnew_vec: s-only density path (numat=1 for a clean check) ---
    {
        int save_numat = numat;
        numat = 1;
        p.assign(4, 0.0);
        p[1] = 0.5;                       // packed diag of atom 1
        nfirst = {0, 1}; nlast = {0, 1};  // nao = 0 (s-only)
        nijbo.assign(2, std::vector<int>(2, 0));
        std::vector<double> vb(3, 0.0);
        bpnew_vec(vb);
        // cosurf[1]=(0.5,0,0), atom1 at (0,0,0): r=2 -> v1 = -0.5*2
        // cosurf[2]=(2,0,0): r=0.5 -> v2 = -0.5*0.5
        chk_rel(vb[1], -0.5 * 2.0, 1e-12, "bpnew_vec v1 (s path)");
        chk_rel(vb[2], -0.5 * 0.5, 1e-12, "bpnew_vec v2 (s path)");
        numat = save_numat;
        nfirst = {0, 1, 2}; nlast = {0, 1, 2};
        nijbo.assign(3, std::vector<int>(3, 0));
    }

    // --- get_bvec: nao=0 -> w[1]=1/r only ---
    {
        double x1[3] = {0, 0, 0}, x2[3] = {1.4, 0, 0};
        double w[46]; for (int i = 0; i < 46; ++i) w[i] = 0.0;
        get_bvec(x1, x2, 0, 1, w);
        chk_rel(w[1], 1.0 / 1.4, 1e-12, "get_bvec s: w1 = 1/r");
        chk(w[2] == 0.0, "get_bvec s: w2 untouched");
        // sp element: dipole/quadrupole terms on axis
        for (int i = 0; i < 46; ++i) w[i] = 0.0;
        get_bvec(x1, x2, 1, 1, w);
        double r = 1.0 / 1.4;
        double dip = dd[1] * a0 * r * r;
        double quad = (a0 * qq[1]) * (a0 * qq[1]) * r * r * r;
        chk_rel(w[1], r, 1e-12, "get_bvec sp: w1 = 1/r");
        chk_rel(w[2], -1.0 * dip, 1e-12, "get_bvec sp: w2 = -dx*dip (dx=-1)");
        chk_rel(w[3], r + (3.0 - 1.0) * quad, 1e-12, "get_bvec sp: w3 = r+2*quad");
        chk_rel(w[5], 0.0, 1e-12, "get_bvec sp: w5 = 0 (off-axis)");
    }

    // --- mult_triangle_vec: packed 3x3 lower triangle ---
    {
        // T = [[2,0,0],[1,3,0],[4,5,6]] packed row-major: 2,1,3,4,5,6
        double t[7] = {0, 2, 1, 3, 4, 5, 6};
        double x[4] = {0, 1, 2, 3};
        double y[4] = {0, 0, 0, 0};
        mult_triangle_vec(t, x, 3, y);
        // symmetric T_sym = [[2,1,4],[1,3,5],[4,5,6]] packed lower row-major
        // y1 = 2*1 + 1*2 + 4*3 = 16
        // y2 = 1*1 + 3*2 + 5*3 = 22
        // y3 = 4*1 + 5*2 + 6*3 = 32
        chk_rel(y[1], 16.0, 1e-12, "mult_triangle y1");
        chk_rel(y[2], 22.0, 1e-12, "mult_triangle y2");
        chk_rel(y[3], 32.0, 1e-12, "mult_triangle y3");
    }

    // --- some_norm ---
    {
        double vv[4] = {0, -3.0, 1.0, 5.0};
        chk_rel(some_norm(vv, 3), 5.0, 1e-12, "some_norm = max abs");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
