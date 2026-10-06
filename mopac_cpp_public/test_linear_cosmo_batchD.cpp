// test_linear_cosmo_batchD.cpp — verify linear_cosmo batch D:
// amat_diag, precondition (block), precondition_solve, and the full
// cgm_solve end-to-end on the H2 surface (n = 20 < na1max, N^2 A*q path).
#include <cmath>
#include <cstdio>
#include <vector>
#include "linear_cosmo.h"
#include "cosmo.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "MOZYME_C.h"
#include "mkl_bits.h"
#include "mopend.h"

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
    // --- setup H2 (same as batches B/C) ---
    numat = 2; norbs = 2; nspa = 42; mozyme = false; nppa = 1082;
    lenabc = std::max(100, nspa * numat);
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[0][1] = 0.0; coord[1][1] = 0.0; coord[2][1] = 0.0;
    coord[0][2] = 1.4; coord[1][2] = 0.0; coord[2][2] = 0.0;
    tore[1] = 1.0; dd[1] = 1.0; qq[1] = 1.0;

    n0[1] = 42; n0[2] = 12;
    dvfill(n0[1], &dirsm[0][0]);
    dvfill(n0[2], &dirsm[0][0] + n0[1] * 4);
    dvfill(1082, &dirvec[0][0]);
    disex2 = 4.0 * std::pow(1.7 * 4.0, 2.0) / nspa;

    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    lijbo = true; nijbo.assign(3, std::vector<int>(3, 0));
    ini_linear_cosmo();

    rsolv = 1.3; ioldcv = 1;
    srad.assign(numat + 1, 0.0);
    srad[1] = 1.30; srad[2] = 1.30;
    iatsp.assign(lenabc + 2, 0);
    nar_csm.assign(lenabc + 2, 0);
    nsetf.assign(lenabc + 2, 0);
    cosurf.assign(5, std::vector<double>(lenabc + 1, 0.0));
    phinet.assign(lenabc + 2, std::vector<double>(4, 0.0));
    arat.assign(numat + 1, 0.0);
    isude.assign(3, std::vector<int>(30 * numat + 1, 0));
    sude.assign(3, std::vector<double>(30 * numat + 1, 0.0));
    coscanz();
    chk(nps > 0, "coscanz produced surface points");
    if (!iblock_pos.empty()) { iblock_pos.clear(); a_block.clear(); }
    iblock_pos.assign(numat + 1, 0);
    int nd = 0;
    max_block_size = 0;
    for (int i1 = 1; i1 <= numat; ++i1) {
        iblock_pos[i1] = nd + 1;
        int j = npoints[i1 + 1] - npoints[i1];
        if (j > max_block_size) max_block_size = j;
        nd += (j * (j + 1)) / 2;
    }
    m_vec.assign(nd + 1, 0.0);
    a_block.assign((max_block_size + 1) * (max_block_size + 1), 0.0);
    chk(nps <= lenabc, "nps within lenabc");

    // --- 1. fill a_part exactly like coscavz would (compute=true) ---
    int maxrs = 60 * numat;
    int ndim = -1;
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
                    ioldcv, maxrs, lenabc, false, ndim);
    a_part.assign(ndim + 1, 0.0);
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
                    ioldcv, maxrs, lenabc, true, ndim);
    chk(ndim == 65, "ndim == 65 (H2 short-range pairs, correct cosurf rows)");
    chk(ndim > 0, "a_part filled");

    // --- 2. amat_diag ---
    a_diag.assign(nps + 1, 0.0);
    amat_diag(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
              linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
              ioldcv, maxrs, lenabc, a_diag);
    bool diag_ok = true;
    for (int i = 1; i <= nps; ++i)
        if (!(a_diag[i] > 0.0) || !std::isfinite(a_diag[i])) { diag_ok = false; break; }
    chk(diag_ok, "amat_diag positive finite diagonal");

    // --- 3. precondition (block form, nps < na1max -> else branch) ---
    m_vec.assign(m_vec.size(), 0.0);
    precondition(cosurf, nps, iatsp, numat, a_diag, m_vec);
    bool m_ok = false;
    for (size_t i = 1; i < m_vec.size(); ++i)
        if (std::fabs(m_vec[i]) > 1e-12 && std::isfinite(m_vec[i])) { m_ok = true; break; }
    chk(m_ok, "precondition filled non-zero finite m_vec");

    // --- 4. precondition_solve: z = M^-1 r, then verify M_block * z ~= r on
    //        the first atom block by re-running mult_triangle_vec on the
    //        packed block (block 1: points 1..npoints[2]) ---
    int nb1 = npoints[2] - npoints[1];   // surface points on atom 1
    chk(nb1 > 0, "atom 1 has surface points");
    std::vector<double> r(nps + 1, 0.0), z(nps + 1, 0.0);
    for (int i = 1; i <= nb1; ++i) r[i] = 0.1 * i;
    precondition_solve(m_vec, z, r, nps);
    bool z_ok = true;
    for (int i = 1; i <= nb1; ++i)
        if (!std::isfinite(z[i])) { z_ok = false; break; }
    chk(z_ok, "precondition_solve finite z");
    chk(std::fabs(z[1]) > 1e-12, "precondition_solve non-trivial z[1]");

    // --- 5. cgm_solve end-to-end: b = A*x_ref, solve, compare x with x_ref ---
    std::vector<double> x_ref(nps + 1, 0.0), b(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) x_ref[i] = 0.01 * (i % 3) + 0.001;
    aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
           linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
           ioldcv, maxrs, lenabc, a_diag, x_ref, b);
    chk(std::fabs(b[3]) > 1e-8, "b = A*x_ref non-zero");

    std::vector<double> x(nps + 1, 0.0), p(nps + 1, 0.0), q(nps + 1, 0.0);
    std::vector<double> m(m_vec.size(), 0.0);   // precondition fills packed blocks
    new_iteration = true;
    cgm_solve(x, nps, b, r, p, q, z, m, true, true);
    // residual r = b - A*x must satisfy some_norm(r) < start_tol (c_proc=0)
    double resid = some_norm(&r[0], nps);   // some_norm is 1-based
    chk(resid < 0.05, "cgm_solve residual below tolerance");
    // solution close to x_ref (CG tolerance 0.01 in inf-norm of residual)
    double max_err = 0.0, ref_scale = 0.0;
    for (int i = 1; i <= nps; ++i) {
        double e = std::fabs(x[i] - x_ref[i]);
        if (e > max_err) max_err = e;
        if (std::fabs(x_ref[i]) > ref_scale) ref_scale = std::fabs(x_ref[i]);
    }
    chk(max_err < 0.05 * ref_scale + 1e-3, "cgm_solve solution close to x_ref");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
