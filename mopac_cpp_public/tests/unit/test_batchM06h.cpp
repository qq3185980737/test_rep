// test_batchM06h.cpp  M06 deriv (internal-coordinate gradient driver) +
// jcarin (Cartesian/internal Jacobian by finite difference).
// Scene: 2 atoms, 1 internal variable (bond length H1-H2 along x), closed
// shell (halfe=false) -> dcart path; NDDO kernels are stubs so dcart gives
// zero Cartesian derivatives. jcarin is tested numerically: forward step on
// the bond gives exactly (0,0,0,step,0,0) in the Jacobian row, so the
// mapped internal gradient equals dxyz(atom2,x) after the 1/step rescale.
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "cosmo_C.h"
#include "deriv.h"
#include "funcon_C.h"
#include "jcarin.h"
#include "molkst_C.h"
#include "molmec_C.h"
#include "MOZYME_C.h"
#include "mxm.h"
#include "parameters_C.h"

static int n_pass = 0, n_fail = 0;
#define CHECK(c, msg)                                     \
  do {                                                    \
    if (c) {                                              \
      ++n_pass;                                           \
      std::printf("  [PASS] %s\n", msg);                  \
    } else {                                              \
      ++n_fail;                                           \
      std::printf("  [FAIL] %s\n", msg);                  \
    }                                                     \
  } while (0)
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a) - (b)) < (tol), msg)

using namespace common_arrays_C;
using namespace molkst_C;
using namespace parameters_C;

// F90 NDDO kernels not in the 2016 tree: stubs keep the pipeline linkable.
extern "C" void diat_(int*, int*, double*, double* smat) {
    for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
    *kr = 1; w[1] = 0.0; *enuc = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
    for (int i = 0; i < 45; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double*, int*, int*) {
    *enuc = 0.0;
}
// post-SCF PM6 correction kernels not yet ported; unused in this scene
// (DH_correction false) but required at link time.
extern "C" void dftd3_() {}
extern "C" void H_bonds4_() {}
extern "C" void energy_corr_hh_rep_() {}
extern "C" void disp_DnX_() {}
extern "C" void PM6_DH_Dispersion_() {}
extern "C" void PM6_DH_H_bond_corrections_() {}
extern "C" void print_post_scf_corrections_() {}

// PM6 correction kernels not yet ported; unused in this scene (method_pm6/
// method_pm7 false) but required at link time by dcart.
void chrge_for_MOZYME(const double*, double*, int, const int*) {}
double nsp2_atom_correction(const std::vector<std::vector<double>>&, int,
                            int, int, int) {
    return 0.0;
}
double Si_O_H_bond_correction(const std::vector<std::vector<double>>&, int,
                              int, int) {
    return 0.0;
}

static void setup_common() {
    numcal = 1;
    norbs = 2; numat = 2; natoms = 2; mpack = 3; nvar = 1; ndep = 0;
    l1u = 0; l2u = 0; l3u = 0; l123 = 1; id = 0; n2elec = 0;
    nclose = 1; nopen = 1; fract = 2.0; gnorm = 0.0;
    cutofp = 9.0;
    keywrd = " ";
    use_ref_geo = false;
    method_pm6 = false;
    method_pm7 = false;
    N_3_present = false;
    Si_O_H_present = false;
    density = 0.0;
    mozyme = false;
    funcon_C::fpc_9 = 14.4; funcon_C::ev = 27.21; funcon_C::a0 = 0.529177;
    nfirst.assign({0, 1, 2});
    nlast.assign({0, 1, 2});
    nat.assign({0, 1, 1});
    p.assign(10, 0.0);
    pa.assign(10, 0.0);
    pb.assign(10, 0.0);
    tvec.assign(4, std::vector<double>(4, 0.0));
    nbonds.assign(5, 0);
    ibonds.assign(5, std::vector<int>(5, 0));
    geo.assign(4, std::vector<double>(4, 0.0));
    geoa.assign(4, std::vector<double>(4, 0.0));
    coord.assign(4, std::vector<double>(4, 0.0));
    labels.assign(4, 0);
    na.assign(4, 0);
    nb.assign(4, 0);
    nc.assign(4, 0);
    step = 0.0;
    for (int i = 0; i <= 4; ++i) tore[i] = 1.0;
    molmec_C::nnhco = 0;
    cosmo_C::useps = false;
    numat_old = 0;
    loc.assign(3, std::vector<int>(4, 0));
    loc[1][1] = 2; loc[2][1] = 1;  // atom 2, bond length
    // z-matrix: atom 1 at origin, atom 2 at bond length 1.5 along x.
    geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
    geo[1][2] = 1.5; geo[2][2] = 0.0; geo[3][2] = 0.0;
    eigs.assign(3, 0.0);
    eigs[1] = -0.5; eigs[2] = 0.3;
    eigb.assign(3, 0.0);
    dxyz.assign(7, 0.0);
}

static void test_jcarin_numeric() {
    std::printf("Test 1: jcarin bond-length Jacobian + mxm mapping\n");
    const double step = 1e-7;
    double xparam[3] = {0.0, 1.5, 0.0};
    double b[8] = {0.0};
    int ncol = 0;
    jcarin(xparam, step, false, b, ncol, 1, 1);
    CHECK(ncol == 6, "jcarin ncol = 3*numat*l123");
    // Row 1: forward coord minus central point = (0,0,0, step, 0, 0).
    CHK_D(b[3], step, 1e-15, "jcarin d(atom2,x)/dq = step");
    CHK_D(b[0], 0.0, 1e-20, "jcarin atom1 x zero");
    CHK_D(b[1], 0.0, 1e-20, "jcarin atom1 y zero");
    CHK_D(b[5], 0.0, 1e-20, "jcarin atom2 z zero");
    // Map a synthetic Cartesian gradient: only atom2-x nonzero.
    double dxyz6[7] = {0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0};
    double gradnt[2] = {0.0, 0.0};
    mxm(b, 1, dxyz6, ncol, &gradnt[1], 1);
    CHK_D(gradnt[1], step * 2.5, 1e-15, "mxm mapped gradient = step*2.5");
}

static void test_deriv_full() {
    std::printf("Test 2: deriv end-to-end (dcart path, stub kernels)\n");
    std::vector<std::vector<double>> geo3(4, std::vector<double>(4, 0.0));
    geo3[1][1] = 0.0; geo3[2][1] = 0.0; geo3[3][1] = 0.0;
    geo3[1][2] = 1.5; geo3[2][2] = 0.0; geo3[3][2] = 0.0;
    std::vector<double> gradnt(4, 0.0);
    gradnt[1] = 0.0;
    deriv(geo3, gradnt);
    CHECK(std::isfinite(gradnt[1]), "deriv gradnt finite");
    CHK_D(gradnt[1], 0.0, 1e-8, "deriv bond gradient zero (zero dxyz)");
    // Second call with a new numcal reinitializes; must stay stable.
    numcal = 2;
    gradnt[1] = 0.0;
    deriv(geo3, gradnt);
    CHK_D(gradnt[1], 0.0, 1e-8, "deriv rerun stable");
    CHECK(std::isfinite(cosine), "deriv cosine finite");
}

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06h batch - deriv / jcarin\n");
    setup_common();
    test_jcarin_numeric();
    setup_common();
    test_deriv_full();
    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
