// test_batchM04l.cpp  M04 density_for_MOZYME (density from compressed LMOs).
// One LMO over two 2-orbital atoms: cocc {1,0,0,1} (atom1 orb1=1, atom2
// orb2=1). nijbo: (1,1)->0 diag block, (2,1)->4 off-diag block, (2,2)->3.
// Hand-computed P (pre-spin): p[1]=1 (c1^2), p[6]=1 (c4^2), p[7]=1 (c1*c4).
// mode=1 doubles (spinfa=2); mode=-1 partpin subtraction then -2;
// mode=2 half partpin then 1.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "density_for_MOZYME.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

// Bare globals (MOZYME_C.cpp, outside namespace).
extern std::vector<int> ncf, ncocc, nncf, icocc;
extern std::vector<double> cocc;

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
  std::printf("M04l batch - density_for_MOZYME\n");
  molkst_C::numat = 2;
  molkst_C::keywrd = " ";
  MOZYME_C::lijbo = true;
  MOZYME_C::iorbs.assign({0, 2, 2});
  ncocc.assign({0, 0});
  ncf.assign({0, 2});
  nncf.assign({0, 0});
  icocc.assign({0, 1, 2});
  cocc.assign(5, 0.0);
  cocc[1] = 1.0;
  cocc[4] = 1.0;
  MOZYME_C::nijbo.assign(3, std::vector<int>(3, -1));
  MOZYME_C::nijbo[1][1] = 0; MOZYME_C::nijbo[2][1] = 4; MOZYME_C::nijbo[2][2] = 3;
  MOZYME_C::nijbo[1][2] = -1;

  // mode=1: full density, spinfa=2.
  std::vector<double> p(9, 0.0);
  std::vector<double> partpin(9, 0.0);
  density_for_MOZYME(p, 1, 1, partpin);
  CHK_D(p[1], 2.0, 1e-12, "dens p(1)=2 (c1^2 * 2)");
  CHK_D(p[6], 2.0, 1e-12, "dens p(6)=2 (c4^2 * 2)");
  CHK_D(p[7], 2.0, 1e-12, "dens p(7)=2 (c1*c4 * 2)");
  CHK_D(p[3], 0.0, 1e-12, "dens p(3)=0 (empty diag)");

  // mode=0: same full density.
  for (auto& x : p) x = 0.0;
  density_for_MOZYME(p, 0, 1, partpin);
  CHK_D(p[1], 2.0, 1e-12, "dens mode0 p(1)=2");

  // mode=2: partial in partpin (0.5 scale), add rest, spinfa=1.
  partpin[1] = 1.0; partpin[6] = 1.0; partpin[7] = 1.0;
  for (auto& x : p) x = 0.0;
  density_for_MOZYME(p, 2, 1, partpin);
  CHK_D(p[1], 1.5, 1e-12, "dens mode2 p(1)=0.5+1");
  CHK_D(p[7], 1.5, 1e-12, "dens mode2 p(7)=0.5+1");

  // mode=-1: subtract partpin (0.5 scale), spinfa=-2.
  for (auto& x : p) x = 0.0;
  density_for_MOZYME(p, -1, 1, partpin);
  CHK_D(p[1], -1.0, 1e-12, "dens mode-1 p(1)=(1-0.5)*-2");
  CHK_D(p[7], -1.0, 1e-12, "dens mode-1 p(7)=(1-0.5)*-2");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
