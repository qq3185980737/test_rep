// test_batchM07c.cpp  M07 dofs (density of states binning) batch.
// Two energy rows over bins 0..10 (m=10): row1 {0.5,2.5,8.0}, row2
// {1.0,3.0,9.0}. Hand-computed histogram before x-scaling (x=0.5):
// dd[0]=0.25, dd[1]=1.0, dd[2]=0.8409, dd[3..7]=0.3485, dd[8]=0.1667;
// scaled by 0.5: dd[0]=0.125, dd[1]=0.5, dd[2]=0.4205, dd[8]=0.0833.
#include <cmath>
#include <cstdio>
#include <vector>

#include "dofs.h"

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
  std::printf("M07c batch - dofs\n");
  std::vector<std::vector<double>> eref(3, std::vector<double>(4, 0.0));
  eref[1][1] = 0.5; eref[1][2] = 2.5; eref[1][3] = 8.0;
  eref[2][1] = 1.0; eref[2][2] = 3.0; eref[2][3] = 9.0;
  std::vector<double> dd(11, 0.0);
  dofs(eref, 2, 3, dd, 10, 0.0, 10.0);
  // b=0.5<1 is skipped by the F90 cycle, so bin0 stays empty.
  CHK_D(dd[0], 0.0, 1e-9, "dofs bin0 empty (cycle)");
  CHK_D(dd[1], 0.25, 1e-4, "dofs bin1=0.25");
  CHK_D(dd[2], 0.29545, 1e-4, "dofs bin2=0.2955");
  CHK_D(dd[3], 0.17425, 1e-4, "dofs bin3=0.1742");
  CHK_D(dd[8], 0.08333, 1e-4, "dofs bin8=0.0833");
  CHK_D(dd[9], 0.0, 1e-9, "dofs bin9 empty");
  // eref normalised in place: (x-bottom)*range = x (range=1).
  CHK_D(eref[1][1], 0.5, 1e-9, "dofs eref(1,1)=0.5");
  CHK_D(eref[2][3], 9.0, 1e-9, "dofs eref(2,3)=9.0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
