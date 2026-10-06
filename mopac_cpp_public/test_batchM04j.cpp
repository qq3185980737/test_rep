// test_batchM04j.cpp  M04 kab_for_MOZYME (K-block two-center Fock add) batch.
// sum(m) = sum_i pk(i)*w(kkind(i,m)); then ia>ja half-triangle j=(j1(j1-1))/2
// adds -sum(m) at j+j2; ia<ja at (j2(j2-1))/2+j1; ia==ja at f(i)-sum(i)*0.5.
// All-ones inputs give sum(m)=16 for every m. Hand-computed offsets below.
#include <cmath>
#include <cstdio>
#include <vector>

#include "kab_for_MOZYME.h"

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
  std::printf("M04j batch - kab_for_MOZYME\n");
  double pk[17], w[101], f[64];
  for (int i = 1; i <= 16; ++i) pk[i] = 1.0;
  for (int k = 1; k <= 100; ++k) w[k] = 1.0;
  for (int k = 0; k < 64; ++k) f[k] = 0.0;

  // ia>ja: j1=3,4,5,6; j=(3,6,10,15); j2=1..4 -> j3 4,5,6,7 | 7,8,9,10 |
  // 11,12,13,14 | 16,17,18,19. f[7] hit twice (j1=3,j2=4 and j1=4,j2=1).
  kab_for_MOZYME(3, 1, pk, w, f);
  CHK_D(f[4], -16.0, 1e-12, "kab f(4)=-16 (j1=3,j2=1)");
  CHK_D(f[7], -32.0, 1e-12, "kab f(7)=-32 (two hits)");
  CHK_D(f[19], -16.0, 1e-12, "kab f(19)=-16 (last)");
  CHK_D(f[3], 0.0, 1e-15, "kab f(3) untouched (before block)");

  // ia<ja: j1=1..4, j2=3..6 -> m order 4,7,11,16 | 5,8,12,17 | 6,9,13,18 |
  // 7,10,14,19. f[7] again hit twice (m=2 and m=13).
  for (int k = 0; k < 64; ++k) f[k] = 0.0;
  kab_for_MOZYME(1, 3, pk, w, f);
  CHK_D(f[4], -16.0, 1e-12, "kab ia<ja f(4)=-16 (m=1)");
  CHK_D(f[7], -32.0, 1e-12, "kab ia<ja f(7)=-32 (two hits)");
  CHK_D(f[19], -16.0, 1e-12, "kab ia<ja f(19)=-16 (last)");
  CHK_D(f[1], 0.0, 1e-15, "kab ia<ja f(1) untouched");

  // ia==ja: f(i) -= sum(i)*0.5 = -8 for i=1..16.
  for (int k = 0; k < 64; ++k) f[k] = 0.0;
  kab_for_MOZYME(2, 2, pk, w, f);
  CHK_D(f[1], -8.0, 1e-12, "kab diag f(1)=-8 (0.5 factor)");
  CHK_D(f[16], -8.0, 1e-12, "kab diag f(16)=-8");
  CHK_D(f[17], 0.0, 1e-15, "kab diag f(17) untouched");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
