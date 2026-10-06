// test_batchC17.cpp — semantic tests for the C17 batch (lbfgs + big_swap full
// translations).  Verifies:
//  1. setulb drives Rosenbrock to convergence via the full L-BFGS-B FG loop
//  2. big_swap store/extract roundtrip preserves MOZYME arrays
//  3. l_control add/remove keyword editing
//  4. build_active_site SET branch builds the set from loc
//  5. get_pars with no calib.dat leaves nloop = 0
//  6. lbfsav restart save/restore roundtrip
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "big_swap.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "lbfgs.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace molkst_C;
using namespace MOZYME_C;

static int failures = 0;
#define CHECK(cond, msg)                                  \
  do {                                                    \
    if (cond) {                                           \
      std::printf("  [PASS] %s\n", msg);                  \
    } else {                                              \
      std::printf("  [FAIL] %s\n", msg);                  \
      ++failures;                                         \
    }                                                     \
  } while (0)

// ---------------------------------------------------------------------------
// 1. setulb full L-BFGS-B run on Rosenbrock, n=2, m=4.
// ---------------------------------------------------------------------------
static void test_setulb_rosenbrock() {
  std::printf("Test 1: setulb Rosenbrock convergence\n");
  int n = 2, m = 4;
  int niwa = 3 * n;
  int nwa = 2 * n * m + 4 * n + 11 * m * m + 8 * m;
  std::vector<double> x((size_t)n + 1, 0.0), l((size_t)n + 1, 0.0);
  std::vector<double> u((size_t)n + 1, 0.0), g((size_t)n + 1, 0.0);
  std::vector<double> wa((size_t)nwa + 1, 0.0), dsave(45, 0.0);
  int nbd[4] = {0, 0, 0, 0};
  std::vector<int> iwa((size_t)niwa + 1, 0);
  std::vector<int> lsave(5, 0), isave(45, 0);
  x[1] = -1.2;
  x[2] = 1.0;
  l[1] = -1.e10;
  l[2] = -1.e10;
  u[1] = 1.e10;
  u[2] = 1.e10;
  std::string task = "START";
  std::string csave = " Unused";
  double f = 0.0;
  int nfg = 0;
  int guard = 0;
  while (true) {
    setulb(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], 1.e7, 1.e-5, wa, iwa,
           task, -1, csave, lsave, isave, dsave);
    if (task.size() >= 2 && task.substr(0, 2) == "FG") {
      double x1 = x[1], x2 = x[2];
      f = 100.0 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1.0 - x1) * (1.0 - x1);
      g[1] = -400.0 * x1 * (x2 - x1 * x1) - 2.0 * (1.0 - x1);
      g[2] = 200.0 * (x2 - x1 * x1);
      ++nfg;
    } else if (task.substr(0, 11) == "CONVERGENCE") {
      break;
    } else if (task.substr(0, 5) == "ERROR") {
      break;
    } else {
      ++guard;
      if (guard > 2000) break;
    }
  }
  double err = std::fabs(x[1] - 1.0) + std::fabs(x[2] - 1.0);
  std::printf("    x* = (%.6f, %.6f)  f = %.3e  nfg = %d  task = %s\n",
              x[1], x[2], f, nfg, task.c_str());
  CHECK(err < 0.05, "setulb converges to Rosenbrock minimum (1,1)");
  CHECK(f < 1.e-3, "setulb final function value small");
  CHECK(task.substr(0, 11) == "CONVERGENCE", "setulb reports CONVERGENCE");
}

// ---------------------------------------------------------------------------
// 2. big_swap store/extract roundtrip.
// ---------------------------------------------------------------------------
static void test_big_swap_roundtrip() {
  std::printf("Test 2: big_swap store/extract roundtrip\n");
  int saved_numat = numat;
  numat = 2;
  nbonds.assign(3, 0);
  nbonds[1] = 1;
  nbonds[2] = 1;
  ibonds.assign(16, std::vector<int>(3, 0));
  ibonds[1][1] = 2;
  ibonds[1][2] = 1;
  geo.assign(4, std::vector<double>(3, 0.0));
  geo[1][1] = 0.5;
  geo[2][1] = 1.5;
  geo[3][1] = 2.5;
  geo[1][2] = -0.5;
  xparam.assign(5, 0.0);
  xparam[1] = 1.25;
  dxyz.assign(7, 0.0);
  dxyz[2] = -3.5;
  icocc.assign(3, 0);
  icocc[1] = 4;
  icvir.assign(3, 0);
  icvir[1] = 5;
  cocc.assign(3, 0.0);
  cocc[1] = 6.5;
  cvir.assign(3, 0.0);
  cvir[1] = 7.5;
  p.assign(3, 0.0);
  p[1] = 6.5;
  icocc_dim = 2;
  mpack = 33;

  big_swap(0, 1);  // store system 1

  // Corrupt current arrays.
  geo[1][1] = 999.0;
  xparam[1] = -99.0;
  nbonds[1] = 0;
  dxyz[2] = 123.0;
  icocc[1] = 0;
  cocc[1] = -1.0;

  big_swap(1, 1);  // extract system 1

  CHECK(geo[1][1] == 0.5 && geo[2][1] == 1.5, "geo restored");
  CHECK(xparam[1] == 1.25, "xparam restored");
  CHECK(nbonds[1] == 1, "nbonds restored");
  CHECK(dxyz[2] == -3.5, "dxyz restored");
  CHECK(icocc[1] == 4, "icocc restored");
  CHECK(cocc[1] == 6.5, "cocc restored");
  CHECK(std::fabs(pa[1] - 0.5 * 6.5) < 1.e-12, "pa = 0.5*p after extract");
  numat = saved_numat;
}

// ---------------------------------------------------------------------------
// 3. l_control add/remove.
// ---------------------------------------------------------------------------
static void test_l_control() {
  std::printf("Test 3: l_control keyword editing\n");
  keywrd = " GEO_REF=foo GNORM=0.5  LET";
  l_control("GNORM=0.5", 9, -1);
  CHECK(keywrd.find("GNORM") == std::string::npos, "l_control removes GNORM");
  l_control("XYZ", 3, 1);
  CHECK(keywrd.find("XYZ") != std::string::npos, "l_control adds XYZ");
  l_control("XYZ", 3, -1);
  CHECK(keywrd.find("XYZ") == std::string::npos, "l_control removes XYZ again");
  // 2-arg form
  l_control("ABC", 1);
  CHECK(keywrd.find("ABC") != std::string::npos, "l_control 2-arg adds ABC");
}

// ---------------------------------------------------------------------------
// 4. build_active_site SET branch.
// ---------------------------------------------------------------------------
static void test_build_active_site_set() {
  std::printf("Test 4: build_active_site SET branch\n");
  keywrd = " LOCATE-TS(SET2)";
  int saved_nvar = nvar;
  nvar = 2;
  loc.assign(3, std::vector<int>(10, 0));
  loc[1][1] = 3;
  loc[1][2] = 7;
  std::vector<int> active_site(201, 0);
  std::vector<int> ninsite(4, 0);
  build_active_site(active_site, ninsite);
  CHECK(ninsite[1] == 2, "SET branch ninsite(1) == 2");
  CHECK(active_site[1] == 3 && active_site[2] == 7, "SET branch atoms from loc");
  nvar = saved_nvar;
}

// ---------------------------------------------------------------------------
// 5. get_pars with no calib.dat.
// ---------------------------------------------------------------------------
static void test_get_pars_nofile() {
  std::printf("Test 5: get_pars without calib.dat\n");
  std::remove("calib.dat");
  std::vector<double> stresses(21, 0.0), gradients(21, 0.0);
  std::vector<double> relscf(21, 0.0), cutoff(21, 0.0);
  int nloop = 99;
  get_pars(stresses, gradients, relscf, cutoff, nloop);
  CHECK(nloop == 0, "get_pars leaves nloop = 0 without input file");
}

// ---------------------------------------------------------------------------
// 6. lbfsav restart save/restore roundtrip.
// ---------------------------------------------------------------------------
static void test_lbfsav_roundtrip() {
  std::printf("Test 6: lbfsav restart roundtrip\n");
  restart_fn = "test_restart.bin";
  std::vector<double> wa(6, 0.0);
  std::vector<int> iwa(6, 0);
  wa[1] = 1.5;
  wa[5] = -2.25;
  iwa[1] = 3;
  iwa[4] = -7;
  std::string task = "START";
  std::string csave = " Unused";
  std::vector<int> lsave(5, 0), isave(45, 0);
  std::vector<double> dsave(30, 0.0);
  lsave[1] = 1;
  isave[3] = 11;
  dsave[4] = 2.5;
  int nstep = 7;
  double escf = -123.45;
  double tt0 = 4321.5;

  lbfsav(tt0, 1, wa, 5, iwa, 4, task, csave, lsave, isave, dsave, nstep, escf);

  // Clear all state.
  wa.assign(6, 0.0);
  iwa.assign(6, 0);
  lsave.assign(5, 0);
  isave.assign(45, 0);
  dsave.assign(30, 0.0);
  nstep = 0;
  escf = 0.0;
  std::string task2 = " Unused";

  int saved_moperr = moperr ? 1 : 0;
  moperr = false;
  lbfsav(0.0, 0, wa, 5, iwa, 4, task2, csave, lsave, isave, dsave, nstep, escf);
  CHECK(!moperr, "lbfsav restore succeeds");
  CHECK(std::fabs(wa[1] - 1.5) < 1.e-12 && std::fabs(wa[5] - (-2.25)) < 1.e-12,
        "lbfsav wa restored");
  CHECK(iwa[1] == 3 && iwa[4] == -7, "lbfsav iwa restored");
  CHECK(nstep == 7, "lbfsav nstep restored");
  CHECK(std::fabs(escf - (-123.45)) < 1.e-9, "lbfsav escf restored");
  CHECK(std::fabs(dsave[4] - 2.5) < 1.e-12, "lbfsav dsave restored");
  moperr = (saved_moperr != 0);
  std::remove("test_restart.bin");
}

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);
  std::printf("==== batch C17 tests: lbfgs.cpp + big_swap.cpp ====\n");
  test_setulb_rosenbrock();
  test_big_swap_roundtrip();
  test_l_control();
  test_build_active_site_set();
  test_get_pars_nofile();
  test_lbfsav_roundtrip();
  if (failures == 0) {
    std::printf("ALL PASS\n");
    return 0;
  }
  std::printf("%d FAILURE(S)\n", failures);
  return 1;
}