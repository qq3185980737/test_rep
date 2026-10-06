// test_linear_cosmo_batchF.cpp — verify linear_cosmo batch F:
// addnucz (nuclear charges into qdenet) and addfckz (Fock correction +
// solvation energy, N^2 branch).  H2 surface, setup same as batches B/C/D/E.
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

int main() {
    // --- setup H2 (same as batches B/C/D/E) ---
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
    int maxrs = 60 * numat;

// --- A-part + preconditioner (as in batch D steps 1-3) ---
    int ndim = -1;
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
                    false, ndim);
a_part.assign(ndim + 1, 0.0);
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc,
                    true, ndim);
    a_diag.assign(nps + 1, 0.0);
    amat_diag(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
              linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm, ioldcv, maxrs, lenabc, a_diag);
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
precondition(cosurf, nps, iatsp, numat, a_diag, m_vec);

    // --- arrays for addfckz / addnucz ---
    fepsi = 78.4;
    p.assign(64, 0.0);
    for (int m = 1; m <= 45; ++m) p[m] = 0.01 * m;
    f.assign(16, 0.0);
    lm61 = 9;
    ipiden.assign(lm61 + 1, 0);
    for (int i = 1; i <= lm61; ++i) ipiden[i] = i;
    gden.assign(lm61 + 1, 0.0);
    for (int i = 1; i <= lm61; ++i) gden[i] = 0.5;
    idenat.assign(numat + 1, 0);
    idenat[1] = 1; idenat[2] = 2;
    qscnet.assign(lenabc + 2, std::vector<double>(4, 0.0));
    qdenet.assign(lm61 + 1, std::vector<double>(4, 0.0));
    qscat.assign(numat + 1, 0.0);
    r_vec.assign(nps + 1, 0.0); p_vec.assign(nps + 1, 0.0);
    q_vec.assign(nps + 1, 0.0); z_vec.assign(nps + 1, 0.0);

// ---- 1. addnucz ----
    new_surface = true;
addnucz(phinet, qscnet, qdenet);
    bool nuc_ok = true;
    for (int i = 1; i <= nps; ++i)
        if (phinet[i][1] != 0.0 || qscnet[i][1] != 0.0) { nuc_ok = false; break; }
    for (int i = 1; i <= lm61; ++i) {
        double want = 0.0;
        for (int a = 1; a <= numat; ++a)
            if (idenat[a] == i) want = tore[nat[a]];
        if (std::fabs(qdenet[i][1] - want) > 1e-12) { nuc_ok = false; break; }
    }
    chk(nuc_ok, "addnucz zeros channels and loads nuclear charges");

// ---- 2. addfckz (new_surface=true) ----
addfckz();
    chk(new_surface == false, "addfckz clears new_surface");
    chk(std::fabs(phinet[1][2]) > 1e-12 && std::isfinite(phinet[1][2]),
        "addfckz potential phi(2) non-zero finite");
    bool sol_ok = true;
    for (int i = 1; i <= nps; ++i) {
        if (!std::isfinite(qscnet[i][2]) || std::fabs(qscnet[i][2]) < 1e-12) { sol_ok = false; break; }
        if (std::fabs(qscnet[i][3] - (qscnet[i][1] + qscnet[i][2])) > 1e-12) { sol_ok = false; break; }
    }
    chk(sol_ok, "addfckz induced charges q(2)/q(3) consistent");
    double qsum = 0.0;
    for (int i = 1; i <= nps; ++i) qsum += qscnet[i][3];
    chk(std::fabs((qscat[1] + qscat[2]) - qsum) < 1e-9 * (1.0 + std::fabs(qsum)),
        "addfckz qscat sums the atomic charges");
    chk(std::isfinite(ediel) && std::isfinite(solv_energy), "addfckz ediel/solv_energy finite");
    chk(std::fabs(solv_energy) > 1e-12, "addfckz solv_energy non-zero");
    chk(std::fabs(f[1]) > 1e-12, "addfckz updated Fock f(1) (S element)");
    bool qd_ok = true;
    for (int i = 1; i <= lm61; ++i) {
        if (std::fabs(qdenet[i][2] - gden[i] * p[ipiden[i]]) > 1e-12) { qd_ok = false; break; }
        if (std::fabs(qdenet[i][3] - (qdenet[i][2] + qdenet[i][1])) > 1e-12) { qd_ok = false; break; }
    }
    chk(qd_ok, "addfckz qdenet from density matrix");

    // ---- 3. addfckz again (new_surface=false warm path) ----
    double q2_before = qscnet[1][2];
    addfckz();
    chk(std::fabs(qscnet[1][2]) > 1e-12 && std::isfinite(qscnet[1][2]),
        "addfckz warm path keeps finite charges");
    (void)q2_before;

    std::fprintf(stderr, "ALL %d CHECKS PASS (%d FAIL)\n", g_checks, g_fail);
    (void)0;
    return g_fail ? 1 : 0;
}
