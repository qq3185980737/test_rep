// test_batchM07a.cpp  M07 linmin (line minimisation) batch.
// compfg stub: f(x) = (x1-1)^2 + x2^2, gradient 2(x1-1), 2*x2.
// Search direction pvect = (1,0) so x1 = alpha.
// Scene 1 (numcal=1, keywrd=" "): Thiel alpha gives alpha<2 with okf
// -> goto 190; hand-computed funct=0.36 (f at x=0.4), xparam={0.4,0},
// alpha=0.4, ic=1.
// Scene 2 (numcal=2, keywrd=" NOTH"): delta1=0.1 so full parabolic loop
// runs and converges to x1=1, funct->0.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "linmin.h"
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
#define CHK_D(a, b, tol, msg)                                            \
  do {                                                                   \
    double va = (a), vb = (b);                                           \
    if (std::fabs(va - vb) <= (tol)) {                                   \
      ++n_pass;                                                          \
      std::printf("  [PASS] %s (%.6g)\n", msg, va);                      \
    } else {                                                             \
      ++n_fail;                                                          \
      std::printf("  [FAIL] %s: got %.6g want %.6g\n", msg, va, vb);     \
    }                                                                    \
  } while (0)

void compfg(const std::vector<double>& xp, bool, double& f, bool,
            std::vector<double>& g, bool) {
  f = (xp[1] - 1.0) * (xp[1] - 1.0) + xp[2] * xp[2];
  g[1] = 2.0 * (xp[1] - 1.0);
  g[2] = 2.0 * xp[2];
}

int main() {
  std::printf("M07a batch - linmin\n");
  double xparam[3] = {0.0, 0.0, 0.0};
  double pvect[3] = {0.0, 1.0, 0.0};
  pvect[0] = 1.0; pvect[1] = 0.0;

  molkst_C::numcal = 1;
  molkst_C::keywrd = " ";
  double alpha = 1.0, funct = 1.0;
  bool okf = false;
  int ic = 0;
  linmin(xparam, alpha, pvect, 2, funct, okf, ic, 1.0);
  CHK_D(funct, 0.36, 1e-9, "linmin early-exit funct=0.36");
  CHK_D(xparam[0], 0.4, 1e-9, "linmin early-exit x(1)=0.4");
  CHK_D(alpha, 0.4, 1e-9, "linmin early-exit alpha=0.4");
  CHECK(ic == 1, "linmin ic=1 (energy unchanged)");
  CHECK(okf, "linmin okf=true");

  // Scene 2: NOTH -> delta1=0.1, full parabolic loop converges to x=1.
  molkst_C::numcal = 2;
  molkst_C::keywrd = " NOTH";
  xparam[0] = xparam[1] = xparam[2] = 0.0;
  pvect[0] = 1.0; pvect[1] = 0.0; pvect[2] = 0.0;
  alpha = 1.0; funct = 1.0; okf = false; ic = 0;
  linmin(xparam, alpha, pvect, 2, funct, okf, ic, 1.0);
  CHECK(funct < 1e-3, "linmin converged funct<1e-3");
  CHK_D(xparam[0], 1.0, 1e-2, "linmin converged x(1)=1");
  CHECK(okf, "linmin conv okf=true");
  CHECK(alpha > 0.0, "linmin alpha>0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
