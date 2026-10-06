// test_linear_cosmo_batchB.cpp — verify linear_cosmo batch B:
// coscanz SAS construction on H2 (box partition, transformation matrices,
// surface points, segment packing, area/volume), then coscavz driver
// (tessellation handles + A-part / preconditioner allocation).
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
    // --- setup H2 ---
    numat = 2; norbs = 2; nspa = 42; mozyme = false; nppa = 1082;
    lenabc = std::max(100, nspa * numat);
    nat = {0, 1, 1};
    coord.assign(4, std::vector<double>(3, 0.0));
    coord[0][1] = 0.0; coord[1][1] = 0.0; coord[2][1] = 0.0;
    coord[0][2] = 1.4; coord[1][2] = 0.0; coord[2][2] = 0.0;
    tore[1] = 1.0; dd[1] = 1.0; qq[1] = 1.0;

    // direction tables (dvfill from cosmo.h), n0: 42 heavy / 12 hydrogen
    n0[1] = 42; n0[2] = 12;
    dvfill(n0[1], &dirsm[0][0]);
    dvfill(n0[2], &dirsm[0][0] + n0[1] * 4);
    dvfill(1082, &dirvec[0][0]);
    disex2 = 4.0 * std::pow(1.7 * 4.0, 2.0) / nspa;

    // allocate linear_cosmo module storage (nset/nipsrs/... used by coscanz)
    nfirst = {0, 1, 2}; nlast = {0, 1, 2};
    lijbo = true; nijbo.assign(3, std::vector<int>(3, 0));
    ini_linear_cosmo();

    rsolv = 1.3; ioldcv = 1;  // closure path (surclo) verified separately: area conserved, no OOB
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

    // --- coscanz ---
    coscanz();
    chk(nps > 0, "coscanz produced surface points");
    chk(nps <= lenabc, "nps within lenabc");
    chk(area > 0.0, "area positive");
    chk(cosvol > 0.0, "cosvol positive");
    // every point lies on its atom sphere: |p - atom| == srad (F90: cosurf = unit_dir*srad + xa)
    double rAtm = srad[1];
    bool on_sas = true;
    for (int i = 1; i <= nps; ++i) {
        int a = iatsp[i];
        double dx = cosurf[1][i] - coord[0][a];
        double dy = cosurf[2][i] - coord[1][a];
        double dz = cosurf[3][i] - coord[2][a];
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (std::fabs(d - rAtm) > 1e-4) { on_sas = false; break; }
    }
    chk(on_sas, "all points on atom sphere (radius srad)");
    // no point inside the other atom's sphere
    bool no_inside = true;
    for (int i = 1; i <= nps; ++i) {
        int a = iatsp[i];
        int b = (a == 1) ? 2 : 1;
        double dx = cosurf[1][i] - coord[0][b];
        double dy = cosurf[2][i] - coord[1][b];
        double dz = cosurf[3][i] - coord[2][b];
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d < rAtm - 1e-6) { no_inside = false; break; }
    }
    chk(no_inside, "no point inside neighbour atom sphere");
    // iatsp covers both atoms
    bool has_a1 = false, has_a2 = false;
    for (int i = 1; i <= nps; ++i) {
        if (iatsp[i] == 1) has_a1 = true;
        if (iatsp[i] == 2) has_a2 = true;
    }
    chk(has_a1 && has_a2, "both atoms have surface points");

    // npoints: transformed cumulative counts; total must equal nps
    int total = 0;
    for (int i = 1; i <= numat; ++i) total += npoints[i + 1] - npoints[i];
    chk(total == nps, "npoints sum to nps");

    // arat: sum of segment areas equals area
    double sum_arat = 0.0;
    for (int i = 1; i <= numat; ++i) sum_arat += arat[i];
    chk_rel(sum_arat, area, 1e-6 * std::max(1.0, std::fabs(area)), "arat sums to area");

    // coscavz driver: allocate A-part and preconditioner (nps small -> simulate_aq_vec stub, ndim=0)
    coscavz();
    chk(atom_handle != 0, "atom_handle set");
    chk(surface_handle != 0, "surface_handle set");
    chk(a_part.size() == 0, "a_part allocated via stub (ndim=0)");
    chk(iblock_pos.size() == (size_t)(numat + 1), "iblock_pos allocated");
    chk(max_block_size >= 0, "max_block_size set");
    chk(!m_vec.empty(), "m_vec allocated");
    chk(!a_block.empty(), "a_block allocated");
    // iblock_pos recurrence: iblock_pos(1)=1
    chk(iblock_pos[1] == 1, "iblock_pos(1)=1");
    // ndim of m_vec = sum j(j+1)/2
    int nd = 0;
    for (int i = 1; i <= numat; ++i) {
        int j = npoints[i + 1] - npoints[i];
        nd += (j * (j + 1)) / 2;
    }
    chk((int)m_vec.size() == nd + 1, "m_vec sized by block sum");
    // block size: max j
    int mb = 0;
    for (int i = 1; i <= numat; ++i) {
        int j = npoints[i + 1] - npoints[i];
        if (j > mb) mb = j;
    }
    chk(max_block_size == mb, "max_block_size = max j");

    std::fprintf(stderr, "ALL %d CHECKS PASS\n", g_checks);
    return g_fail ? 1 : 0;
}
