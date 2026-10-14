// test_batchM07h.cpp  M07 dfpsav (DFP restart dump/restore) batch.
// Scenario A: mdfp(9)=1 dump -> file header norbs/numat written,
//             geout/den_in_out(1) called, no mopend.
// Scenario B: mdfp(9)=0 restore with numcal bumped -> xparam/funct1/gd/
//             hesinv recovered, den_in_out(0) called.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "dfpsav.h"
#include "ef_C.h"
#include "maps_C.h"
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

static int g_geout = 0, g_mopend = 0, g_prtgra = 0, g_den = 0;
extern "C" void geout(int) { ++g_geout; }
void mopend(const char*) { ++g_mopend; }
void prtgra() { ++g_prtgra; }
void den_in_out(int mode) { ++g_den; }

int main() {
  std::printf("M07h batch - dfpsav\n");
  using namespace common_arrays_C;
  using namespace molkst_C;
  using namespace maps_C;

  chanel_C::restart_fn = "dfpsav_test.bin";
  std::remove("dfpsav_test.bin");
  numcal = 1;
  nvar = 2;
  norbs = 4;
  numat = 1;
  natoms = 1;
  latom = 0;
  keywrd = " DFPSAV-TEST";
  koment = "comment";
  title = "title";

  geo.resize(4);
  for (int j = 1; j <= 3; ++j) geo[j].assign(2, 0.0);
  geoa.resize(4);
  for (int j = 1; j <= 3; ++j) geoa[j].assign(2, 0.0);
  na.assign(2, 0);

  std::vector<double> xparam(3, 0.0), gd(3, 0.0), xlast(3, 0.0);
  xparam[1] = 1.5; xparam[2] = 2.5;
  gd[1] = 0.5; gd[2] = 0.7;
  xlast[1] = 1.4; xlast[2] = 2.4;
  grad.assign(3, 0.0);
  grad[1] = 0.2; grad[2] = 0.3;
  hesinv.assign(4, 0.0);
  hesinv[1] = 0.11; hesinv[2] = 0.22; hesinv[3] = 0.33;
  double totime = 12.5, funct1 = -12.5;
  std::vector<int> mdfp(10, 0);
  std::vector<double> xdfp(10, 0.0);
  xdfp[1] = 1.0; xdfp[9] = 9.0;

  // Scenario A: dump
  mdfp[9] = 1;
  dfpsav(totime, xparam, gd, xlast, funct1, mdfp, xdfp);
  CHECK(g_geout == 1, "dump calls geout once (mdfp(9)==1)");
  CHECK(g_den == 1, "dump calls den_in_out(1)");
  CHECK(g_mopend == 0, "dump no mopend");
  CHECK(std::ifstream("dfpsav_test.bin", std::ios::binary).good(),
        "restart file created");
  {
    std::ifstream f("dfpsav_test.bin", std::ios::binary);
    int h1 = 0, h2 = 0;
    f.read(reinterpret_cast<char*>(&h1), sizeof(int));
    f.read(reinterpret_cast<char*>(&h2), sizeof(int));
    CHECK(h1 == 4 && h2 == 1, "file header norbs=4 numat=1");
  }

  // Scenario B: restore (bump numcal so first==true)
  numcal = 2;
  std::vector<double> r_xp(3, -99.0), r_gd(3, -99.0), r_xl(3, -99.0);
  double r_totime = -1.0, r_funct = -1.0;
  std::vector<double> r_xdfp(10, -1.0);
  gd.assign(3, 0.0);
  grad.assign(3, 0.0);
  hesinv.assign(4, 0.0);
  mdfp.assign(10, 0);
  mdfp[9] = 0;
  dfpsav(r_totime, r_xp, r_gd, r_xl, r_funct, mdfp, r_xdfp);
  CHECK(g_den == 2, "restore calls den_in_out(0)");
  CHECK(g_mopend == 0, "restore no mopend (valid file)");
  CHECK(std::abs(r_xp[1] - 1.5) < 1e-12, "xparam(1) restored 1.5");
  CHECK(std::abs(r_xp[2] - 2.5) < 1e-12, "xparam(2) restored 2.5");
  CHECK(std::abs(r_gd[2] - 0.7) < 1e-12, "gd(2) restored 0.7");
  CHECK(std::abs(r_funct + 12.5) < 1e-12, "funct1 restored -12.5");
  CHECK(std::abs(r_totime - 12.5) < 1e-12, "totime restored 12.5");
  CHECK(std::abs(hesinv[3] - 0.33) < 1e-12, "hesinv(3) restored 0.33");
  CHECK(std::abs(r_xdfp[9] - 9.0) < 1e-12, "xdfp(9) restored 9.0");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
