// test_linear_cosmo_batchG.cpp — verify linear_cosmo batch G:
// fock_dir_int / fock_mult_int / fock_far_int / sphere_f_multipoles /
// am1dft_solve.  H2 surface setup shared with batches B..F.
#include <cmath>
#include <complex>
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

// Independent recomputation helpers (share only afmm_mod::get_legendre).
static double ref_local_exp(double x0, double y0, double z0,
                            const afmm_mod::CArr& psi, const afmm_mod::RArr& y_norm,
                            int p, double cx, double cy, double cz, bool far) {
    afmm_mod::RArr pmn;
    double dx = cx - x0, dy = cy - y0, dz = cz - z0;
    double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    double cos_t = dz / d, phi = std::atan2(dy, dx);
    afmm_mod::get_legendre(p, cos_t, pmn);
    double t = 0.0;
    for (int j = 0; j <= p; ++j) {
        std::complex<double> tc(0.0, 0.0);
        for (int k = 1; k <= j; ++k) {
            double mul = far ? 1.0 / std::pow(d, j + 1) : std::pow(d, j);
            tc += y_norm(k, j) * pmn(k, j) * mul * psi(k, j) *
                  std::exp(std::complex<double>(0.0, k * phi));
        }
        double mul0 = far ? 1.0 / std::pow(d, j + 1) : std::pow(d, j);
        t += 2.0 * std::real(tc) + pmn(0, j) * mul0 * std::real(psi(0, j));
    }
    return t;
}

int main() {
    // --- setup H2 (same as batches B..F) ---
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
    qscnet.assign(lenabc + 2, std::vector<double>(4, 0.0));
    qdenet.assign(10, std::vector<double>(4, 0.0));
    arat.assign(numat + 1, 0.0);
    isude.assign(3, std::vector<int>(30 * numat + 1, 0));
    sude.assign(3, std::vector<double>(30 * numat + 1, 0.0));
    coscanz();
    chk(nps > 0, "coscanz produced surface points");
    int maxrs = 60 * numat;

    // --- A-part + preconditioner (for am1dft_solve) ---
    int ndim = -1;
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
                    ioldcv, maxrs, lenabc, false, ndim);
    a_part.assign(ndim + 1, 0.0);
    simulate_aq_vec(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
                    linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
                    ioldcv, maxrs, lenabc, true, ndim);
    a_diag.assign(nps + 1, 0.0);
    amat_diag(coord, srad, numat, cosurf, nps, nar_csm, nsetf,
              linear_cosmo::nset, rsc, nipsrs, iatsp, linear_cosmo::tm,
              ioldcv, maxrs, lenabc, a_diag);
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

    // --- batch G scratch ---
    fepsi = 78.4;
    p.assign(64, 0.0);
    for (int m = 1; m <= 45; ++m) p[m] = 0.01 * m;
    iatom_pos.assign(numat + 1, 0);
    iatom_pos[1] = 0; iatom_pos[2] = 0;
    lm61 = 9;
    ipiden.assign(lm61 + 1, 0);
    for (int i = 1; i <= lm61; ++i) ipiden[i] = i;
    gden.assign(lm61 + 1, 0.0);
    for (int i = 1; i <= lm61; ++i) gden[i] = 0.5;
    idenat.assign(numat + 1, 0);
    idenat[1] = 1; idenat[2] = 2;
    qscat.assign(numat + 1, 0.0);
    for (int i = 1; i <= nps; ++i) qscnet[i][2] = 0.01 * i;

    // column-major atom coordinates for the AFMM callbacks
    std::vector<double> cflat(3 * numat, 0.0);
    for (int i = 1; i <= numat; ++i) {
        cflat[(i - 1) * 3 + 0] = coord[0][i];
        cflat[(i - 1) * 3 + 1] = coord[1][i];
        cflat[(i - 1) * 3 + 2] = coord[2][i];
    }
    int inda[3] = {0, 1, 2};   // 1-based ind: ind[1]=1, ind[2]=2

    // ---- 1. fock_dir_int (same-box: atoms 1,2 vs their own points) ----
    {
        std::vector<double> r(8, 0.0);
        fock_dir_int(inda, 2, inda, 2, cflat.data(), 3, cflat.data(), r.data(), true);
        // independent recomputation: v(1) = sum over all surface points of
        // get_bvec(...)w(1) * qscnet(j,2), written to r[0+1]
        double ref = 0.0;
        for (int i = 1; i <= numat; ++i) {
            double xa[3] = {coord[0][i], coord[1][i], coord[2][i]};
            double w[46];
            for (int j = 1; j <= nps; ++j) {
                double xp[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
                get_bvec(xp, xa, 0, 1, w);
                ref += w[1] * qscnet[j][2];
            }
        }
        chk_rel(r[1], ref, 1e-9 * (1.0 + std::fabs(ref)), "fock_dir_int same r(1)");
        bool rest_zero = true;
        for (int m = 2; m < 8; ++m) if (r[m] != 0.0) { rest_zero = false; break; }
        chk(rest_zero, "fock_dir_int writes only ii+1");
    }
    // ---- 1b. fock_dir_int (!same: bidirectional two boxes, ind1={1}, ind2={2}) ----
    {
        int i1[2] = {0, 1}, i2[2] = {0, 2};
        std::vector<double> r(8, 0.0);
        fock_dir_int(i1, 1, i2, 1, cflat.data(), 3, cflat.data(), r.data(), false);
        double ref = 0.0;
        // (1) atom 1 with points of atom 2 (npoints[2]..npoints[3]-1)
        double xa[3] = {coord[0][1], coord[1][1], coord[2][1]};
        double w[46];
        for (int j = npoints[2]; j <= npoints[3] - 1; ++j) {
            double xp[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
            get_bvec(xp, xa, 0, 1, w);
            ref += w[1] * qscnet[j][2];
        }
        // (2) atom 2 with points of atom 1 (npoints[1]..npoints[2]-1)
        xa[0] = coord[0][2]; xa[1] = coord[1][2]; xa[2] = coord[2][2];
        for (int j = npoints[1]; j <= npoints[2] - 1; ++j) {
            double xp[3] = {cosurf[1][j], cosurf[2][j], cosurf[3][j]};
            get_bvec(xp, xa, 0, 1, w);
            ref += w[1] * qscnet[j][2];
        }
        chk_rel(r[1], ref, 1e-9 * (1.0 + std::fabs(ref)), "fock_dir_int !same bidirectional");
    }

    // ---- 2. fock_mult_int (local expansion) ----
    {
        afmm_mod::CArr psi;
        afmm_mod::RArr y_norm;
        psi(0, 0) = std::complex<double>(0.3, 0.0);
        for (int j = 1; j <= 3; ++j) {
            psi(0, j) = std::complex<double>(0.1 * j, 0.0);
            for (int k = 1; k <= j; ++k) {
                psi(k, j) = std::complex<double>(0.05 * k, 0.02 * j);
                y_norm(k, j) = 0.1 * k + 0.01 * j;
            }
        }
        const double x0 = 5.0, y0 = 0.0, z0 = 1.0;
        std::vector<double> r(8, 0.0);
        afmm_mod::RArr pmn_g;
        fock_mult_int(inda, 2, psi, 3, x0, y0, z0, y_norm, pmn_g,
                      cflat.data(), 3, cflat.data(), r.data());
        double ref = 0.0;
        for (int i = 1; i <= numat; ++i)
            ref += ref_local_exp(x0, y0, z0, psi, y_norm, 3,
                                 coord[0][i], coord[1][i], coord[2][i], false);
        chk_rel(r[1], ref, 1e-9 * (1.0 + std::fabs(ref)), "fock_mult_int r(1)");
    }

    // ---- 3. fock_far_int (multipole expansion) ----
    {
        afmm_mod::CArr psi;
        afmm_mod::RArr y_norm;
        psi(0, 0) = std::complex<double>(0.3, 0.0);
        for (int j = 1; j <= 3; ++j) {
            psi(0, j) = std::complex<double>(0.1 * j, 0.0);
            for (int k = 1; k <= j; ++k) {
                psi(k, j) = std::complex<double>(0.05 * k, 0.02 * j);
                y_norm(k, j) = 0.1 * k + 0.01 * j;
            }
        }
        const double x0 = 5.0, y0 = 0.0, z0 = 1.0;
        std::vector<double> r(8, 0.0);
        afmm_mod::RArr pmn_g;
        fock_far_int(inda, 2, psi, 3, x0, y0, z0, y_norm, pmn_g,
                     cflat.data(), 3, cflat.data(), r.data());
        double ref = 0.0;
        for (int i = 1; i <= numat; ++i)
            ref += ref_local_exp(x0, y0, z0, psi, y_norm, 3,
                                 coord[0][i], coord[1][i], coord[2][i], true);
        chk_rel(r[1], ref, 1e-9 * (1.0 + std::fabs(ref)), "fock_far_int r(1)");
    }

    // ---- 4. sphere_f_multipoles ----
    {
        afmm_mod::CArr psi;   // zero-init
        afmm_mod::RArr y_norm, pmn;
        for (int j = 1; j <= 3; ++j)
            for (int k = 0; k <= j; ++k) y_norm(k, j) = 0.05 * k + 0.01 * j;
        const double x0 = 3.0, y0 = 0.5, z0 = -0.5;
        sphere_f_multipoles(x0, y0, z0, inda, 2, pmn, y_norm,
                            cflat.data(), 3, cflat.data(), psi, 3);
        // psi(0,0) = sum of qq over all surface points
        double qsum = 0.0;
        for (int j = 1; j <= nps; ++j) qsum += qscnet[j][2];
        chk_rel(std::real(psi(0, 0)), qsum, 1e-12, "sphere_f_multipoles psi(0,0) charge");
        chk(std::fabs(std::imag(psi(0, 0))) < 1e-12, "sphere_f_multipoles psi(0,0) real");
        // independent recomputation of psi(m1,n1) and conjugation symmetry
        bool sym_ok = true;
        for (int n1 = 1; n1 <= 3; ++n1)
            for (int m1 = 1; m1 <= n1; ++m1)
                if (std::abs(psi(-m1, n1) - std::conj(psi(m1, n1))) > 1e-12) sym_ok = false;
        chk(sym_ok, "sphere_f_multipoles conjugation symmetry");
        // ref psi(2,2) element: sum over points of exp(-2i phi)*y_norm(-2,2)*pmn(-2,2)*r^2*qq
        std::complex<double> ref22(0.0, 0.0);
        for (int j = 1; j <= nps; ++j) {
            double dx = cosurf[1][j] - x0, dy = cosurf[2][j] - y0, dz = cosurf[3][j] - z0;
            double rr = std::sqrt(dx * dx + dy * dy + dz * dz);
            double cos_t = dz / rr, phi = std::atan2(dy, dx);
            afmm_mod::RArr pm;
            afmm_mod::get_legendre(3, cos_t, pm);
            std::complex<double> tc = std::exp(std::complex<double>(0.0, -2.0 * phi));
            ref22 += tc * (y_norm(-2, 2) * pm(-2, 2) * rr * rr * qscnet[j][2]);
        }
        chk_rel(std::abs(psi(2, 2) - ref22), 0.0, 1e-9,
                "sphere_f_multipoles psi(2,2) element");
    }

    // ---- 5. am1dft_solve (cgm_solve wrapper, warm path) ----
    {
        std::vector<double> qsct(nps + 1, 0.0);
        std::vector<double> phit(nps + 1, 0.0);
        for (int i = 1; i <= nps; ++i) phit[i] = 0.5 + 0.01 * i;
        am1dft_solve(qsct, phit, nps);
        bool finite_nonzero = true;
        for (int i = 1; i <= nps; ++i)
            if (!std::isfinite(qsct[i]) || std::fabs(qsct[i]) < 1e-12) { finite_nonzero = false; break; }
        chk(finite_nonzero, "am1dft_solve produces finite non-zero charges");
    }

    std::fprintf(stderr, "ALL %d CHECKS PASS (%d FAIL)\n", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
