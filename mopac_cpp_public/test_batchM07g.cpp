// test_batchM07g.cpp  M07 drc (DRC/IRC driver) batch.
// Single atom, keywrd " CYCLES=5" -> IRC path (half=0), compfg stub
// f=x1^2/grad=2x1, prtdrc stub increments itemp_1 (Fortran jloop).
// Driver must run cycles until itemp_1==5 then exit restoring iw0.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "drc.h"
#include "elemts_C.h"
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

static int g_compfg = 0, g_prtdrc = 0;
static double g_t = 0.0;

double second(int n) {
  if (n == 1) return 0.0;
  g_t += 0.001;
  return g_t;
}
double reada(const std::string& s, int istart) {
  // " CYCLES=5" -> parse digits after position.
  std::size_t p = s.find("CYCLES=");
  if (p == std::string::npos) return 0.0;
  p += 7;
  return std::atof(s.c_str() + p);
}
void compfg(const std::vector<double>& xp, bool, double& f, bool,
            std::vector<double>& g, bool) {
  ++g_compfg;
  f = xp[1] * xp[1];
  g[1] = 2.0 * xp[1];
  g[2] = 0.0;
  g[3] = 0.0;
}
void prtdrc(double, const std::vector<double>&,
            const std::vector<std::vector<double>>&, double, double, double,
            const std::vector<double>&, const std::vector<std::vector<int>>&,
            int, bool) {
  ++g_prtdrc;
  ++molkst_C::itemp_1;
}
void gmetry(std::vector<std::vector<double>>& geo,
            std::vector<std::vector<double>>& coord) {
  // Test stub: deliver the current internal coordinates as cartesian.
  coord = geo;
}
void l_control(const char*, int, int) {}
void mopend(const char*) {}
void den_in_out(int) {}
double ddot(int n, const double* x, int, const double* y, int) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[i] * y[i];
  return s;
}
extern "C" void to_screen_(const char*) {}

int main() {
  std::printf("M07g batch - drc\n");
  using namespace molkst_C;
  using namespace common_arrays_C;
  using namespace elemts_C;

  numcal = 1;
  keywrd = " CYCLES=3";
  natoms = 1;
  numat = 1;
  nvar = 0;
  nopen = 1; nclose = 1;
  moperr = false;
  ndep = 0;
  itemp_1 = 0;
  chanel_C::iw0 = -1;
  tleft = 1e6;
  chanel_C::restart_fn = "drc_restart_test.bin";
  atmass.resize(2);
  atmass[1] = 12.0;
  geo.resize(4);
  for (int j = 1; j <= 3; ++j) geo[j].assign(2, 0.0);
  geo[1][1] = 1.0;
  coord.resize(4);
  for (int j = 1; j <= 3; ++j) coord[j].assign(2, 0.0);
  nat.resize(2, 0);
  nat[1] = 6;
  elemnt.resize(20);
  elemnt[6] = " C ";
  loc.resize(3);
  loc[1].assign(4, 0);
  loc[2].assign(4, 0);

  std::vector<double> startv(3 * 1, 0.0), startk(3 * 1, 0.0);
  drc(startv, startk);

  CHECK(g_compfg == 3, "drc ran 3 compfg calls");
  CHECK(g_prtdrc == 3, "drc ran 3 prtdrc calls");
  CHECK(itemp_1 == 3, "drc itemp_1 reaches 3");
  CHECK(nvar == 3, "drc nvar=3*numat");
  CHECK(chanel_C::iw0 == -1, "drc iw0 restored");
  CHECK(std::abs(escf - 1.0) < 0.01, "drc escf final ~1.0");
  CHECK(gnorm > 0.0, "drc gnorm positive (grad non-zero)");
  CHECK(xparam[1] > 0.99 && xparam[1] < 1.01,
        "drc xparam(1) stays near 1 (tiny steps)");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
