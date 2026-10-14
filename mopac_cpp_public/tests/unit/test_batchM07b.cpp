// test_batchM07b.cpp  M07 locmin (line search on |grad|^2) batch.
// compfg stub: f=(x1-1)^2+x2^2, grad = (2(x1-1), 2x2). Direction (1,0).
// locmin minimises ddot(efs,efs)=4((x1-1)^2+x2^2) along p; minimum at
// x1=1 where gradient vanishes -> ssq->0, x1->1, ncount>=2.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "locmin.h"
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
double ddot(int n, const double* x, int, const double* y, int) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[i] * y[i];
  return s;
}

int main() {
  std::printf("M07b batch - locmin\n");
  molkst_C::numcal = 1;
  molkst_C::keywrd = " ";
  double xparam[3] = {0.0, 0.0, 0.0};
  double p[3] = {0.0, 1.0, 0.0};
  p[0] = 1.0; p[1] = 0.0;
  double efs[3] = {0.0, 0.0, 0.0};
  double ssq = 4.0;   // |grad(0,0)|^2 = 4
  double alf = 1.0;
  int ncount = 0;
  locmin(2, xparam, 2, p, ssq, alf, efs, ncount);
  CHECK(ssq < 1e-4, "locmin ssq->0 (gradient vanished)");
  CHK_D(xparam[0], 1.0, 1e-2, "locmin x1->1");
  CHECK(alf > 0.0, "locmin alf>0");
  CHECK(ncount >= 2, "locmin ncount>=2");
  CHK_D(efs[0], 0.0, 1e-2, "locmin efs restored (grad ~0)");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
