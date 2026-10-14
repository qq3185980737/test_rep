// test_batchM07m.cpp  M07m: ef (P-RFO/QA optimizer, ef.F90) batch.
// Scenario A: nvar=1, x1=1.5 -> one EF step (diag Hessian 1000) moves to
//             x1=1.499, then rmx=0.998 < tol2=1.0 -> iflepo=15.
// Scenario B: nvar=2 at minimum -> gnorm=0 < tol2 -> iflepo=2.
// Scenario C: " RESTART" with missing restart file -> efsav fails -> moperr.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "ef.h"
#include "ef_C.h"
#include "maps_C.h"
#include "molkst_C.h"

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

static int g_compfg = 0, g_mopend = 0;

void compfg(const std::vector<double>& xp, bool, double& f, bool,
            std::vector<double>& g, bool lgrad) {
  ++g_compfg;
  double x1 = xp.size() > 1 ? xp[1] : 0.0;
  double x2 = xp.size() > 2 ? xp[2] : 0.0;
  f = (x1 - 1.0) * (x1 - 1.0) + (x2 - 2.0) * (x2 - 2.0);
  if (lgrad) {
    g[1] = 2.0 * (x1 - 1.0);
    if (g.size() > 2) g[2] = 2.0 * (x2 - 2.0);
  }
}
double ddot(int n, const double* x, int, const double* y, int) {
  double s = 0.0;
  for (int i = 1; i <= n; ++i) s += x[i] * y[i];
  return s;
}
double second(int n) {
  static double t = 0.0;
  if (n == 1) return 0.0;
  t += 0.001;
  return t;
}
double reada(const std::string& s, int istart) {
  std::size_t p = s.find("CYCLES=");
  if (p == std::string::npos) return 0.0;
  p += 7;
  return std::atof(s.c_str() + p);
}
void prttim(double, double& tprt, char& txt) {
  tprt = 0.5;
  txt = 'W';
}
void symtry() {}
void geout(int) {}
void mopend(const char*) {
  ++g_mopend;
  molkst_C::moperr = true;  // Fortran mopend sets the error flag
}
void den_in_out(int) {}
void prtgra() {}
void write_cell(int) {}
void to_screen(const std::string&) {}

int main() {
  std::printf("M07m batch - ef\n");
  std::fflush(stdout);
  using namespace molkst_C;
  using namespace common_arrays_C;
  chanel_C::restart_fn = "ef_restart.bin";
  std::remove("ef_restart.bin");
  numcal = 1;
  tleft = 1e6;
  time0 = 0.0;
  tdump = 100000;
  moperr = false;
  iflepo = 0;
  last = 0;
  nscf = 1;
  natoms = 1;
  numat = 1;
  ndep = 0;
  id = 0;
  norbs = 1;
  nalpha = 0;
  keywrd = "  TEST";
  loc.assign(3, std::vector<int>(4, 0));
  geo.assign(4, std::vector<double>(2, 0.0));
  nc.assign(4, 0);  // ef best-geometry bookkeeping reads nc[1..natoms]

  // --- Scenario A: nvar=1, x1=1.5 ---
  nvar = 1;
  loc[1][1] = 1;  // atom 1
  loc[2][1] = 1;  // coordinate component 1 (x)
  grad.assign(2, 0.0);
  hesinv.clear();
  g_compfg = 0;
  g_mopend = 0;
  std::vector<double> xp(4, 0.0);
  xp[1] = 1.5;
  double f = 99.0;
  std::printf("  [A] calling ef nvar=1\n");
  std::fflush(stdout);
  ef(xp, f);
  std::printf("  [A] returned iflepo=%d f=%.6f x1=%.6f\n", iflepo, f, xp[1]);
  std::fflush(stdout);
  CHECK(iflepo == 15, "A: nvar=1 converges -> iflepo=15");
  CHECK(std::abs(xp[1] - 1.499) < 1e-6,
        "A: x1 -> 1.499 after one EF step");
  // compfg: f = (x1-1)^2 + (x2-2)^2 with x2 fixed at 0 -> (0.499)^2 + 4
  CHECK(std::abs(f - 4.249001) < 1e-6,
        "A: f=(1.499-1)^2+(0-2)^2=4.249001");
  CHECK(g_mopend == 0, "A: no mopend");

  // --- Scenario B: nvar=2 at minimum -> gnorm=0 -> iflepo=2 ---
  nvar = 2;
  loc[1][1] = 1; loc[2][1] = 1;
  loc[1][2] = 1; loc[2][2] = 2;
  numcal = 2;
  moperr = false;
  g_compfg = 0;
  g_mopend = 0;
  grad.assign(3, 0.0);
  hesinv.clear();
  std::vector<double> xp2(4, 0.0);
  xp2[1] = 1.0;
  xp2[2] = 2.0;
  double f2 = 99.0;
  std::printf("  [B] calling ef nvar=2\n");
  std::fflush(stdout);
  ef(xp2, f2);
  std::printf("  [B] returned iflepo=%d f=%.6f\n", iflepo, f2);
  std::fflush(stdout);
  CHECK(iflepo == 2, "B: at minimum -> iflepo=2 (gradient cutoff)");
  CHECK(std::abs(f2) < 1e-9, "B: f=0 at minimum");
  CHECK(g_mopend == 0, "B: no mopend");

  // --- Scenario C: RESTART with missing file -> moperr ---
  nvar = 2;
  numcal = 3;
  moperr = false;
  keywrd = " RESTART";
  g_mopend = 0;
  hesinv.clear();
  std::vector<double> xp3(4, 0.0);
  xp3[1] = 0.5;
  xp3[2] = 1.0;
  double f3 = 99.0;
  ef(xp3, f3);
  CHECK(moperr == true, "C: RESTART missing file -> moperr (efsav fails)");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
