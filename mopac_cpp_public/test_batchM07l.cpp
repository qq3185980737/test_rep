// test_batchM07l.cpp  M07 flepo (BFGS/DFP optimiser) batch.
// Scenario A: nvar=1 -> linmin then nvar==1 exit (iflepo=14, grad(1)=0).
// Scenario B: nvar=2 at minimum -> gnorm<tolerg -> iflepo=2.
// Scenario C: " RESTAR" with missing restart file -> dfpsav read fails
//             -> moperr set, flepo returns early.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "flepo.h"
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

static int g_compfg = 0, g_mopend = 0, g_geout = 0, g_dfpsav = 0;

void compfg(const std::vector<double>& xp, bool, double& f, bool,
            std::vector<double>& g, bool lgrad) {
  ++g_compfg;
  // f = (x1-1)^2 + (x2-2)^2  (nvar-agnostic: missing vars stay 0)
  double x1 = xp.size() > 1 ? xp[1] : 0.0;
  double x2 = xp.size() > 2 ? xp[2] : 0.0;
  f = (x1 - 1.0) * (x1 - 1.0) + (x2 - 2.0) * (x2 - 2.0);
  if (lgrad) {  // gradient written only when requested (Fortran lgrad)
    g[1] = 2.0 * (x1 - 1.0);
    if (g.size() > 2) g[2] = 2.0 * (x2 - 2.0);
  }
}
void linmin(double* xp, double& alpha, double* pv, int nvar, double& funct,
            bool& okf, int& ic, double) {
  // stub: advance alpha*pv, caller recomputes funct via compfg
  for (int i = 1; i <= nvar; ++i) xp[i] += alpha * pv[i];
  funct = 0.0;
  okf = true;
  ic = 2;
}
void dcopy(int n, const double* x, int, double* y, int) {
  for (int i = 1; i <= n; ++i) y[i] = x[i];
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
extern "C" void geout(int) { ++g_geout; }
void mopend(const char*) {
  ++g_mopend;
  molkst_C::moperr = true;  // Fortran mopend sets the error flag
}
void prtgra() {}
void den_in_out(int) {}

int main() {
  std::printf("M07l batch - flepo\n");
  std::fflush(stdout);
  using namespace molkst_C;
  using namespace common_arrays_C;
  chanel_C::restart_fn = "flepo_restart.bin";
  std::remove("flepo_restart.bin");
  numcal = 1;
  grad.assign(3, 0.0);
  hesinv.assign(4, 0.0);
  tleft = 1e6;
  time0 = 0.0;
  tdump = 100000;
  moperr = false;
  iflepo = 0;
  last = 0;
  nscf = 1;
  keywrd = "  TEST";

  // --- Scenario A: nvar=1 ---
  g_compfg = 0;
  std::vector<double> xp(2, 0.0);
  xp[1] = 1.5;
  double f = 99.0;
  std::printf("  [A] calling flepo nvar=1\n"); std::fflush(stdout);
  flepo(xp, 1, f);
  std::printf("  [A] returned iflepo=%d f=%.6f\n", iflepo, f); std::fflush(stdout);
  CHECK(iflepo == 14, "nvar=1: iflepo=14 (single-variable exit)");
  CHECK(std::abs(grad[1]) < 1e-12,
        "nvar=1: grad(1)=0 after exit (compfg lgrad=false)");
  // linmin stub moves x1: 1.5 -> 1.51 (alpha=1.0, pv=0.01); x2 fixed at 0
  // -> f = (1.51-1)^2 + (0-2)^2 = 4.2601  (Fortran semantics)
  CHECK(std::abs(f - 4.2601) < 1e-9,
        "nvar=1: funct recomputed at moved x1=1.51 (f=4.2601)");
  CHECK(g_compfg >= 1, "nvar=1: compfg called (final)");

  // --- Scenario B: nvar=2 at minimum -> gnorm<tolerg -> iflepo=2 ---
  numcal = 2;
  moperr = false;
  g_compfg = 0;
  std::vector<double> xp2(3, 0.0);
  xp2[1] = 1.0; xp2[2] = 2.0;
  double f2 = 99.0;
  std::printf("  [B] calling flepo nvar=2\n"); std::fflush(stdout);
  flepo(xp2, 2, f2);
  std::printf("  [B] returned iflepo=%d f=%.6f\n", iflepo, f2);
  std::fflush(stdout);
  CHECK(iflepo == 2, "nvar=2 at min: iflepo=2 (gradient test)");
  CHECK(std::abs(f2) < 1e-9, "nvar=2: funct ~ 0 at minimum");
  CHECK(g_mopend == 0, "nvar=2: no mopend");

  // --- Scenario C: RESTART with missing file -> moperr ---
  numcal = 3;
  moperr = false;
  keywrd = " RESTAR";
  std::vector<double> xp3(3, 0.0);
  xp3[1] = 0.5; xp3[2] = 1.0;
  double f3 = 99.0;
  flepo(xp3, 2, f3);
  CHECK(moperr == true, "RESTART missing file: moperr set (dfpsav fails)");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
