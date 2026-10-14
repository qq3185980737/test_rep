// test_batchM04k.cpp  M04 local_for_MOZYME (LMO psi^4 localisation) batch.
// Two 2-orbital atoms, two LMOs sharing both atoms: LMO1={1,0,0,1},
// LMO2={0,1,1,0} (1-based c). Hand-computed: aij=2, bij=0 -> ca=sqrt(0.5),
// sa=sqrt(0.5); after one rotation c[1]=c[4]=sqrt(2), c[5..8]=0; second
// sweep xiijj=0 -> totij=0 -> loop exits.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "local_for_MOZYME.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

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

extern void mopend(const std::string& msg) { (void)msg; }
extern void memory_error(const std::string& msg) { (void)msg; }

int main() {
  std::printf("M04k batch - local_for_MOZYME\n");
  molkst_C::norbs = 4;
  molkst_C::numat = 2;
  molkst_C::natoms = 2;
  molkst_C::nelecs = 4;
  MOZYME_C::iorbs.assign({0, 2, 2});

  cocc_dim = 8;
  cocc.assign(9, 0.0);
  cocc[1] = 1.0; cocc[4] = 1.0; cocc[6] = 1.0; cocc[7] = 1.0;
  // LMO1 atoms {1,2}; LMO2 atoms {2,1}; LMO2 starts at c index 5.
  icocc_dim = 4;
  icocc.assign({0, 1, 2, 2, 1});
  ncf.assign({0, 2, 2});
  ncocc.assign({0, 0, 4});
  nncf.assign({0, 0, 2});

  local_for_MOZYME("OCCUPIED");
  CHK_D(cocc[1], 1.4142135623730951, 1e-9, "local c(1)=sqrt(2) after rotation");
  CHK_D(cocc[4], 1.4142135623730951, 1e-9, "local c(4)=sqrt(2)");
  CHK_D(cocc[6], 0.0, 1e-12, "local c(6)=0 (LMO2 atom2 orb1 cleared)");
  CHK_D(cocc[7], 0.0, 1e-12, "local c(7)=0 (LMO2 atom1 orb1 cleared)");
  CHK_D(cocc[2], 0.0, 1e-12, "local c(2) untouched 0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
