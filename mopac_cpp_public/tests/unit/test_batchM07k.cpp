// test_batchM07k.cpp  M07 pinout (LMO density file write/read) batch.
// Scenario A: write (mode=1) with nocc=1/nvir=3 -> density file created,
//             prtlmo not called (no PINOUT keyword).
// Scenario B: read (mode=0) -> ncf/nce/icocc/icvir/iorbs/nbonds/ibonds/
//             cocc/cvir + compressed nncf/nnce/ncocc/ncvir recovered.
// Scenario C: missing file -> mopend.
// Scenario D: PINOUT keyword -> prtlmo called after read.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include "molkst_C.h"
#include "pinout.h"

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

static int g_mopend = 0, g_prtlmo = 0, g_screen = 0;
void mopend(const std::string&) { ++g_mopend; }
void prtlmo() { ++g_prtlmo; }
void to_screen(const std::string&) { ++g_screen; }

// MOZYME module data (bare globals).
extern std::vector<int> nce, ncf, ncvir, nncf, nnce, ncocc, icocc, icvir;
extern std::vector<double> cocc, cvir;
extern int cocc_dim, cvir_dim, icocc_dim, icvir_dim;

int main() {
  std::printf("M07k batch - pinout\n");
  using namespace molkst_C;
  using namespace common_arrays_C;
  chanel_C::density_fn = "pinout_test.den";
  std::remove("pinout_test.den");

  nelecs = 2;            // nocc = 1
  norbs = 4;             // nvir = 3
  numat = 2;
  keywrd = "  TEST";

  // --- write-mode setup (nocc=1, nvir=3) ---
  ncf.assign(3, 0); nce.assign(5, 0);
  nncf.assign(3, 0); nnce.assign(5, 0);
  ncocc.assign(3, 0); ncvir.assign(5, 0);
  icocc.assign(4, 0); icvir.assign(6, 0);
  cocc.assign(6, 0.0); cvir.assign(8, 0.0);
  ncf[1] = 2; nce[1] = 1; nce[2] = 1; nce[3] = 1;
  nncf[1] = 0; nnce[1] = 0; nnce[2] = 1; nnce[3] = 2;
  ncocc[1] = 0; ncvir[1] = 0; ncvir[2] = 2; ncvir[3] = 3;
  icocc[1] = 1; icocc[2] = 2;
  icvir[1] = 1; icvir[2] = 2; icvir[3] = 1;
  MOZYME_C::iorbs.assign(3, 0);
  MOZYME_C::iorbs[1] = 2; MOZYME_C::iorbs[2] = 1;
  cocc[1] = 0.1; cocc[2] = 0.2; cocc[3] = 0.3;   // LMO1: iorbs1+iorbs2=3
  cvir[1] = 0.5; cvir[2] = 0.6;                  // LMO1(vir): iorbs[1]=2
  cvir[3] = 0.7;                                  // LMO2: iorbs[2]=1
  cvir[4] = 0.8; cvir[5] = 0.9;                  // LMO3: iorbs[1]=2
  nbonds.assign(3, 0);
  nbonds[1] = 2; nbonds[2] = 1;
  ibonds.resize(10);
  for (int j = 1; j <= 9; ++j) ibonds[j].assign(3, 0);
  ibonds[1][1] = 2; ibonds[2][1] = 3; ibonds[1][2] = 5;
  icocc_dim = 10; icvir_dim = 10; cocc_dim = 10; cvir_dim = 10;

  // Scenario A: write
  pinout(1);
  CHECK(std::ifstream("pinout_test.den", std::ios::binary).good(),
        "density file created");
  CHECK(g_mopend == 0 && g_prtlmo == 0, "write: no mopend, no prtlmo");

  // Scenario B: read (clear first)
  ncf.assign(3, 0); nce.assign(5, 0);
  nncf.assign(3, -1); nnce.assign(5, -1);
  ncocc.assign(3, -1); ncvir.assign(5, -1);
  icocc.assign(4, 0); icvir.assign(6, 0);
  cocc.assign(6, 0.0); cvir.assign(8, 0.0);
  MOZYME_C::iorbs.assign(3, 0);
  nbonds.assign(3, 0);
  for (int j = 1; j <= 9; ++j) ibonds[j].assign(3, 0);
  keywrd = "  TEST";  // no PINOUT -> no prtlmo
  pinout(0);
  CHECK(ncf[1] == 2 && nce[1] == 1 && nce[2] == 1 && nce[3] == 1,
        "ncf/nce recovered");
  CHECK(icocc[1] == 1 && icocc[2] == 2, "icocc recovered");
  CHECK(icvir[1] == 1 && icvir[2] == 2 && icvir[3] == 1,
        "icvir recovered");
  CHECK(MOZYME_C::iorbs[1] == 2 && MOZYME_C::iorbs[2] == 1,
        "iorbs recovered");
  CHECK(nbonds[1] == 2 && nbonds[2] == 1, "nbonds recovered");
  CHECK(ibonds[1][1] == 2 && ibonds[2][1] == 3 && ibonds[1][2] == 5,
        "ibonds recovered");
  CHECK(nncf[1] == 0 && nnce[1] == 0 && nnce[2] == 1 && nnce[3] == 2,
        "compressed nncf/nnce rebuilt");
  CHECK(std::abs(cocc[1] - 0.1) < 1e-12 && std::abs(cocc[3] - 0.3) < 1e-12,
        "cocc LMO coefficients recovered");
  CHECK(ncocc[1] == 0, "ncocc rebuilt 0");
  CHECK(std::abs(cvir[1] - 0.5) < 1e-12 && std::abs(cvir[5] - 0.9) < 1e-12,
        "cvir LMO coefficients recovered");
  CHECK(ncvir[1] == 0 && ncvir[2] == 2 && ncvir[3] == 3,
        "ncvir rebuilt 0/2/3");
  CHECK(g_mopend == 0 && g_prtlmo == 0, "read: no mopend, no prtlmo");

  // Scenario C: missing file
  std::remove("pinout_test.den");
  int g0 = g_mopend;
  pinout(0);
  CHECK(g_mopend == g0 + 1, "missing file -> mopend once");

  // Scenario D: PINOUT keyword -> prtlmo (rewrite file first)
  keywrd = " PINOUT";
  pinout(1);   // rewrite density file
  pinout(0);   // read -> prtlmo triggered
  CHECK(g_prtlmo >= 1, "PINOUT keyword -> prtlmo called");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
