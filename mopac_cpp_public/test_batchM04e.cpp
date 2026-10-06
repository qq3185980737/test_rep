// test_batchM04e.cpp  M04 rotlmo (rotate p AOs in all LMOs) batch.
// One occupied LMO with a p shell (atom 1) and one virtual LMO with a p
// shell are rotated by a 90-degree matrix; the s part (ka+0) is untouched.
// All expectations hand-computed from sum = cocc(ka+k)*rotvec(k,j).
#include <cmath>
#include <cstdio>
#include <vector>

#include "molkst_C.h"
#include "MOZYME_C.h"
#include "rotlmo.h"

using namespace molkst_C;
using namespace MOZYME_C;

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
  std::printf("M04e batch - rotlmo\n");
  nelecs = 2; norbs = 5;   // nocc=1, nvir=4
  iorbs.assign({0, 4, 2, 2, 2, 2});

  // Occupied LMO 1: atom 1 (4 orbitals), starts at cocc[1].
  ncocc.assign({0, 0, 0});
  nncf.assign({0, 0, 0});
  ncf.assign({0, 1, 0});
  icocc.assign({0, 1, 0});
  cocc.assign(20, 0.0);
  cocc[1] = 0.9;   // s
  cocc[2] = 1.0;   // px
  cocc[3] = 0.0;   // py
  cocc[4] = 0.0;   // pz

  // Virtual LMO 1: atom 1 again, starts at cvir[6].
  ncvir.assign({0, 5, 0, 0, 0, 0});
  nnce.assign({0, 0, 0, 0, 0, 0});
  nce.assign({0, 1, 0, 0, 0, 0});
  icvir.assign({0, 1, 0, 0, 0, 0});
  cvir.assign(20, 0.0);
  cvir[6] = 0.3;   // s
  cvir[7] = 0.0;   // px
  cvir[8] = 0.5;   // py
  cvir[9] = 0.0;   // pz

  // 90-degree rotation about z: x->y, y->-x (rotvec(k,j)).
  // rotvec = [[0,-1,0],[1,0,0],[0,0,1]] row-major.
  double rotvec[9] = {0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  rotlmo(rotvec);

  // Occupied: vec(j) = sum_k cocc(ka+k)*rotvec(k,j).
  // j=1: 1*0 + 0*1 + 0*0 = 0 ; j=2: 1*(-1) = -1 ; j=3: 0.
  CHK_D(cocc[1], 0.9, 1e-12, "rotlmo occ s untouched");
  CHK_D(cocc[2], 0.0, 1e-12, "rotlmo occ px -> 0");
  CHK_D(cocc[3], -1.0, 1e-12, "rotlmo occ py -> -1");
  CHK_D(cocc[4], 0.0, 1e-12, "rotlmo occ pz -> 0");

  // Virtual: (px,py,pz)=(0,0.5,0); vec = v^T*R = 0.5*(row1) = (0.5,0,0).
  CHK_D(cvir[6], 0.3, 1e-12, "rotlmo vir s untouched");
  CHK_D(cvir[7], 0.5, 1e-12, "rotlmo vir px -> +0.5");
  CHK_D(cvir[8], 0.0, 1e-12, "rotlmo vir py -> 0");
  CHK_D(cvir[9], 0.0, 1e-12, "rotlmo vir pz -> 0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
