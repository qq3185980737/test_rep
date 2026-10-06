// test_batchM07e.cpp  M07 reorth (MOZYME LMO re-orthogonalisation) batch.
// Small system: 2 atoms, 2 orbitals each, 1 virtual LMO {0.5,0.5} on atom1,
// 1 occupied LMO {0.6,0.8} on atom1. reorth computes the virt-occ overlap
// sum=0.7 and rotates cvir via adjvec: cvecb(m) -= beta*cveca(m).
// Expected: sumtot=0.7, cvir[1]=0.5-0.7*0.6=0.08, cvir[2]=0.5-0.7*0.8=-0.06;
// ws ends holding the occupied LMO coefficients {0.6,0.8}.
#include <cmath>
#include <cstdio>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "reorth.h"

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

// Bare globals from MOZYME_C.cpp (outside namespace).
extern std::vector<int> nce, ncf, ncocc, ncvir, nnce, nncf, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int nvirtual, noccupied;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;

int main() {
  std::printf("M07e batch - reorth\n");
  molkst_C::numat = 2;
  common_arrays_C::nfirst.resize(3);
  common_arrays_C::nfirst[1] = 1;
  common_arrays_C::nfirst[2] = 3;
  MOZYME_C::iorbs.resize(3);
  MOZYME_C::iorbs[1] = 2;
  MOZYME_C::iorbs[2] = 2;
  MOZYME_C::thresh = 1e-5;

  nvirtual = 1; noccupied = 1;
  cvir_dim = 2; cocc_dim = 2;
  nce.resize(2, 0); nnce.resize(2, 0);
  ncf.resize(2, 0); nncf.resize(2, 0);
  ncocc.resize(2, 0); ncvir.resize(2, 0);
  cvir.resize(3, 0.0); cocc.resize(3, 0.0);
  icvir.resize(2, 0); icocc.resize(2, 0);
  nnce[1] = 0; nce[1] = 1; icvir[1] = 1; ncvir[1] = 0;
  cvir[1] = 0.5; cvir[2] = 0.5;
  nncf[1] = 0; ncf[1] = 1; icocc[1] = 1; ncocc[1] = 0;
  cocc[1] = 0.6; cocc[2] = 0.8;

  double ws[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  reorth(ws);

  CHK_D(cvir[1], 0.08, 1e-9, "reorth cvir(1) rotated to 0.08");
  CHK_D(cvir[2], -0.06, 1e-9, "reorth cvir(2) rotated to -0.06");
  CHK_D(ws[1], 0.6, 1e-9, "reorth ws holds occ LMO (1)=0.6");
  CHK_D(ws[2], 0.8, 1e-9, "reorth ws holds occ LMO (2)=0.8");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
