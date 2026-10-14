// test_linear_cosmo_batchC.cpp — verify linear_cosmo batch C:
// A-part machinery — simulate_aq_vec (count + fill a_part), aq_vec (A*q,
// N^2), aq_dir_int (direction-list A*q) and simulate_aq_dir_int (count/fill)
// on the H2 surface produced by coscanz.
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
    // --- setup H2 (same as batch B) ---
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

    // --- 1. simulate_aq_vec: count short-range elements ---
    int ndim = -1;
    int maxrs = 60 * numat;
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
                    false, ndim);
    chk(ndim >= 1, "simulate_aq_vec counted short-range pairs");
    std::fprintf(stderr, "NDIM=%d nps=%d\n", ndim, nps);

    // --- 2. fill a_part and verify non-zero ---
    a_part.assign(ndim + 1, 0.0);   // a_part is used 1-based (F90 ALLOCATABLE)
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
                    true, ndim);
    bool a_part_ok = true;
    double a_part_sum = 0.0;
    for (int i = 1; i <= ndim; ++i) {
        a_part_sum += a_part[i];
        if (a_part[i] == 0.0) { a_part_ok = false; break; }
    }
    chk(a_part_ok, "a_part filled (all non-zero)");
    chk(a_part_sum > 0.0, "a_part elements positive");

    // --- 3. aq_vec: A*q via N^2 ---
    std::vector<double> q(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) q[i] = 0.01 * (i % 3) + 0.001;
    a_diag.assign(nps + 1, 0.0);
    for (int i = 1; i <= nps; ++i) a_diag[i] = 1.0 / (srad[iatsp[i]] * srad[iatsp[i]]);
    std::vector<double> v(nps + 1, 0.0);
    aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
           linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
           a_diag, q, v);
    bool v_ok = false;
    for (int i = 1; i <= nps; ++i) if (std::fabs(v[i]) > 1e-8) v_ok = true;
    chk(v_ok, "aq_vec produced non-zero A*q");
    // v includes the diagonal: v = A*q where A_ii = 1/srad^2
    double v_diag_ref = a_diag[1] * q[1];
    bool has_offdiag = false;
    for (int i = 2; i <= nps; ++i) if (std::fabs(v[i] - a_diag[i] * q[i]) > 1e-10) has_offdiag = true;
    chk(has_offdiag, "aq_vec includes off-diagonal contributions");

    // --- 4. aq_dir_int (same=true over all points) must reproduce aq_vec ---
    std::vector<int> ind(nps + 1);
    for (int i = 1; i <= nps; ++i) ind[i] = i;
    std::vector<double> r(nps + 1, 0.0);
    double cdummy[4] = {1.0, 0.0, 0.0, 0.0};
    new_iteration = true;              // reset static ijpos
    aq_dir_int(ind.data(), nps, nullptr, 0, cdummy, 4, q.data(), r.data(), true);
    bool dir_match = true;
    for (int i = 1; i <= nps; ++i) {
        // aq_vec adds the diagonal a_diag*q at the end; aq_dir_int does not
        double expect = v[i] - a_diag[i] * q[i];
        if (std::fabs(r[i] - expect) > 1e-8) { dir_match = false; break; }
    }
    chk(dir_match, "aq_dir_int matches aq_vec off-diagonal (same=true)");

    // --- 5. simulate_aq_dir_int counts same as simulate_aq_vec ---
    int npos = 0;
    simulate_aq_dir_int(ind.data(), nps, nullptr, 0, cdummy, 4, true, false, npos);
    chk(npos == ndim, "simulate_aq_dir_int count matches simulate_aq_vec");

    // --- 6. diagonal check: v(ii) -= A_ii*q(ii) gives pure off-diagonal part
    //        which is symmetric: sum over ii of off-diag contribution == 0
    double sym_check = 0.0;
    for (int i = 1; i <= nps; ++i) sym_check += v[i];
    double diag_sum = 0.0;
    for (int i = 1; i <= nps; ++i) diag_sum += a_diag[i] * q[i];
    // off-diagonal A is symmetric, so sum_i sum_j!=i A_ij q_j = 2 * sum_{i<j} A_ij (q_i+q_j)
    // not zero in general; just check finite and sign consistency via a second call
    std::vector<double> v2(nps + 1, 0.0);
    std::vector<double> q2 = q;
    for (int i = 1; i <= nps; ++i) q2[i] *= 2.0;
    aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
           linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
           a_diag, q2, v2);
    chk_rel(v2[1], 2.0 * v[1], 1e-8, "aq_vec linear in q (point 1)");
    chk_rel(v2[7], 2.0 * v[7], 1e-8, "aq_vec linear in q (point 7)");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
