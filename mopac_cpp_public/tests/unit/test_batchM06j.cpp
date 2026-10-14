// test_batchM06j.cpp  M06 compfg (main energy+gradient driver) batch.
// Scenes cover the SCF assembly path, the XFAC/POP population branch, the
// reference-geometry stress term, and the int=false path. Unported kernels
// (cosmo surface, MOZYME SCF, PM6 corrections) are internal stubs in
// compfg.cpp; hcore/iter are stubbed here (SCF numeric closure belongs to the
// M03 subsystem). All assertions are hand-derived from F90 formulas.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "compfg.h"
#include "cosmo_C.h"
#include "derivs_C.h"
#include "funcon_C.h"
#include "molmec_C.h"
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace cosmo_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace molmec_C;
using namespace parameters_C;

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
#define CHK_D(a, b, tol, msg)                                             \
  do {                                                                    \
    double va = (a), vb = (b);                                            \
    if (std::fabs(va - vb) <= (tol)) {                                    \
      ++n_pass;                                                           \
      std::printf("  [PASS] %s (%.9g)\n", msg, va);                       \
    } else {                                                              \
      ++n_fail;                                                           \
      std::printf("  [FAIL] %s: got %.9g want %.9g\n", msg, va, vb);      \
    }                                                                     \
  } while (0)

// NDDO kernels absent from the 2016 tree (extern C stubs, same as M06h).
extern "C" void diat_(int*, int*, double*, double* smat) {
  for (int i = 0; i < 81; ++i) smat[i] = 0.0;
}
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj,
                        double* w, int* kr, double* enuc) {
  for (int i = 0; i < 45; ++i) w[i] = 0.0;
  *enuc = 0.0;
}
extern "C" void elenuc_(int*, int*, int*, int*, double* en) {
  for (int i = 0; i < 45; ++i) en[i] = 0.0;
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double*, int*, int*) {
  for (int i = 0; i < 45; ++i) w[i] = 0.0;
  *enuc = 0.0;
}
// PM6 post-SCF correction family (unported).
extern "C" void dftd3_() {}
extern "C" void H_bonds4_() {}
extern "C" void energy_corr_hh_rep_() {}
extern "C" void disp_DnX_() {}
extern "C" void PM6_DH_Dispersion_() {}
extern "C" void PM6_DH_H_bond_corrections_() {}
extern "C" void print_post_scf_corrections_() {}
// SCF kernels (M03 subsystem gap): control-flow stubs.
void hcore() {}
void iter(double&, bool, bool) {}
void hcore_for_MOZYME() {}
void iter_for_MOZYME(double&) {}

// dcart dependencies (M06c2 pattern).
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
  norbs = 2; nelecs = 2; numat = 2; natoms = 2; mpack = 3;
  nvar = 1; ndep = 0;
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
  lxfac = false;
  moperr = false;
  pressure = 0.0;
  atheat = 0.0;
  enuclr = 0.0;
  elect = 0.0;
  emin = 0.0;
  fpc_9 = 14.4; ev = 27.21; a0 = 0.529177;
  nfirst.assign({0, 1, 2});
  nlast.assign({0, 1, 2});
  nat.assign({0, 1, 1});
  p.assign(10, 0.0);
  pa.assign(10, 0.0);
  pb.assign(10, 0.0);
  pdiag.assign(10, 0.0);
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
  nnhco = 0;
  useps = false;
  numat_old = 0;
  loc.assign(3, std::vector<int>(4, 0));
  loc[1][1] = 2; loc[2][1] = 1;  // atom 2, bond length
  geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
  geo[1][2] = 1.5; geo[2][2] = 0.0; geo[3][2] = 0.0;
  eigs.assign(3, 0.0);
  eigs[1] = -0.5; eigs[2] = 0.3;
  eigb.assign(3, 0.0);
  dxyz.assign(7, 0.0);
}

static void test_compfg_scf_path() {
  std::printf("Test 1: compfg SCF assembly (hcore/iter stubbed)\n");
  setup_common();
  atheat = 0.1; enuclr = 0.567; elect = 0.333;
  std::vector<double> xparam = {0.0, 1.5};
  std::vector<double> grad(2, 0.0);
  double escf = 0.0;
  compfg(xparam, true, escf, false, grad, false);
  // escf = (elect+enuclr)*fpc_9 + atheat = 0.9*14.4 + 0.1
  CHK_D(escf, 0.9 * 14.4 + 0.1, 1e-12, "compfg escf hand value");
  CHK_D(atheat, 0.1, 1e-12, "compfg atheat restored");
  CHK_D(emin, escf, 1e-12, "compfg emin updated");
  CHK_D(geo[1][2], 1.5, 1e-12, "compfg geo from xparam");
}

static void test_compfg_emin_update() {
  std::printf("Test 2: compfg emin keeps lowest energy\n");
  setup_common();
  atheat = 0.1; enuclr = 0.567; elect = 0.333;
  std::vector<double> xparam = {0.0, 1.5};
  std::vector<double> grad(2, 0.0);
  double escf = 0.0;
  compfg(xparam, true, escf, false, grad, false);
  double first = escf;
  elect = 0.1;  // lower electronic energy -> lower escf
  compfg(xparam, true, escf, false, grad, false);
  CHK_D(escf, (0.1 + 0.567) * 14.4 + 0.1, 1e-12, "compfg second escf");
  CHK_D(emin, escf, 1e-12, "compfg emin == new minimum");
  CHECK(escf < first, "compfg escf decreased");
}

static void test_compfg_xfac_pop() {
  std::printf("Test 3: compfg XFAC POP population branch\n");
  setup_common();
  norbs = 0;  // no SCF after population setup
  lxfac = true;
  keywrd = " POP 4 6 10";
  atheat = 0.1; enuclr = 0.5; elect = 0.0;
  std::vector<double> xparam = {0.0, 1.5};
  std::vector<double> grad(2, 0.0);
  double escf = 0.0;
  compfg(xparam, true, escf, false, grad, false);
  // F90: pa = full populations, p = pa, then pa *= 0.5 (density halved).
  CHK_D(pa[1], 2.0, 1e-12, "POP pa(1) = s/2 (F90 halved)");
  CHK_D(pa[3], 1.0, 1e-12, "POP pa(3) = p/6 (F90 halved)");
  CHK_D(pa[15], 1.0, 1e-12, "POP pa(15) = d/10 (F90 halved)");
  CHK_D(p[1], 4.0, 1e-12, "POP p(1) = 2*pa(1)");
  CHK_D(pdiag[1], 4.0, 1e-12, "POP pdiag(1)");
  CHK_D(pdiag[5], 2.0, 1e-12, "POP pdiag(5)");
  CHK_D(escf, 0.5 * 14.4 + 0.1, 1e-12, "POP escf = enuclr*fpc9+atheat");
}

static void test_compfg_stress() {
  std::printf("Test 4: compfg reference-geometry stress term\n");
  setup_common();
  use_ref_geo = true;
  density = 1.0;
  geoa[1][1] = 0.1;   // geo[1][1] = 0 -> d = -0.1
  geoa[1][2] = 1.5;   // match geo[1][2] so only atom-1 x deviates
  atheat = 0.1; enuclr = 0.567; elect = 0.333;
  std::vector<double> xparam = {0.0, 1.5};
  std::vector<double> grad(2, 0.0);
  double escf = 0.0;
  compfg(xparam, true, escf, false, grad, false);
  // stress = sum (geo-geoa)^2 = 0.01; atheat += density*stress = 0.01
  CHK_D(escf, 0.9 * 14.4 + 0.1 + 0.01, 1e-12, "compfg escf with stress");
}

static void test_compfg_int_false() {
  std::printf("Test 5: compfg int=false path (no hcore/iter)\n");
  setup_common();
  atheat = 0.1; enuclr = 0.567; elect = 0.333;
  std::vector<double> xparam = {0.0, 1.5};
  std::vector<double> grad(2, 0.0);
  double escf = 0.0;
  compfg(xparam, false, escf, false, grad, false);
  // int=false: SCF skipped, elect untouched -> same assembly value.
  CHK_D(escf, 0.9 * 14.4 + 0.1, 1e-12, "compfg int=false escf");
}

int main() {
  std::printf("M06j batch - compfg\n");
  test_compfg_scf_path();
  test_compfg_emin_update();
  test_compfg_xfac_pop();
  test_compfg_stress();
  test_compfg_int_false();
  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
