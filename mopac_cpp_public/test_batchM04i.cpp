// test_batchM04i.cpp  M04 jab_for_MOZYME (J-block two-center Fock add) batch.
// Hand-computed from the F90: suma/sumb(i) = sum_j pja/pjb(j)*w(jindx/jjndx)
// then the i5/i6 loop adds sumb(i) to f1 and suma(i) to f2 at the running
// ioff/joff offsets. Scene ia=1, ja=2 verifies the offset ladder.
#include <cmath>
#include <cstdio>
#include <vector>

#include "jab_for_MOZYME.h"

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

int main() {
  std::printf("M04i batch - jab_for_MOZYME\n");
  double pja[17], pjb[17], w[101], f1[64], f2[64];
  for (int j = 1; j <= 16; ++j) {
    pja[j] = 1.0;
    pjb[j] = 1.0;
  }
  for (int k = 1; k <= 100; ++k) w[k] = 1.0;
  for (int k = 0; k < 64; ++k) f1[k] = f2[k] = 0.0;

  // All ones: suma(i)=sumb(i)=16 for every i.
  jab_for_MOZYME(1, 2, pja, pjb, w, f1, f2);
  CHK_D(f1[1], 16.0, 1e-12, "jab f1(1)=16 (i5=1)");
  CHK_D(f1[3], 16.0, 1e-12, "jab f1(3)=16 (i5=2, i6=2)");
  CHK_D(f1[10], 16.0, 1e-12, "jab f1(10)=16 (i5=4, i6=4)");
  CHK_D(f2[3], 16.0, 1e-12, "jab f2(3)=16");
  // Last joff: i5=4, ija=5 -> joff=(5*4)/2+ja-1=11, 4 steps -> 15.
  CHK_D(f2[15], 16.0, 1e-12, "jab f2(15)=16 (last)");
  CHK_D(f1[11], 0.0, 1e-15, "jab f1 beyond block untouched");

  // pjb=0 -> f1 unchanged; pja[j]=j -> suma(i)=136 -> f2(3)=136.
  for (int j = 1; j <= 16; ++j) {
    pja[j] = (double)j;
    pjb[j] = 0.0;
  }
  for (int k = 0; k < 64; ++k) f1[k] = f2[k] = 0.0;
  jab_for_MOZYME(1, 2, pja, pjb, w, f1, f2);
  CHK_D(f2[3], 136.0, 1e-12, "jab f2(3)=sum j = 136");
  CHK_D(f1[1], 0.0, 1e-15, "jab f1 untouched when pjb=0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
