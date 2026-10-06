// test_batchM07j.cpp  M07 forsav (FORCE restart save/restore) batch.
// Scenario A: write (ipt=2) -> den_in_out(1), header time/ipt/refh/numat/
//             norbs, time>1e7 wrapped.
// Scenario B: read (ipt=0) -> all arrays recovered.
// Scenario C: missing file -> mopend.
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "forsav.h"
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

static int g_mopend = 0, g_den = 0;
void mopend(const std::string&) { ++g_mopend; }
void den_in_out(int mode) { ++g_den; }

int main() {
  std::printf("M07j batch - forsav\n");
  using namespace molkst_C;
  chanel_C::restart_fn = "forsav_test.bin";
  std::remove("forsav_test.bin");
  numat = 1;
  norbs = 4;
  const int nvar = 2;

  std::vector<std::vector<double>> deldip(4, std::vector<double>(3, 0.0));
  deldip[1][1] = 0.1; deldip[2][1] = 0.2; deldip[3][1] = 0.3;
  deldip[1][2] = 1.1; deldip[2][2] = 1.2; deldip[3][2] = 1.3;
  std::vector<double> fmatrx(4, 0.0), coord(3, 0.0), evecs(5, 0.0);
  std::vector<double> fconst(3, 0.0);
  fmatrx[1] = 0.5; fmatrx[2] = 0.6; fmatrx[3] = 0.7;
  coord[1] = 1.5; coord[2] = 2.5;
  evecs[1] = 0.9; evecs[2] = 0.8; evecs[3] = 0.7; evecs[4] = 0.6;
  fconst[1] = 3.1; fconst[2] = 3.2;
  double time = 10000005.0, refh = -1.5;
  int ipt = 2, jstart = 7;

  // Scenario A: write (time>1e7 wraps)
  forsav(time, deldip, ipt, fmatrx, coord, nvar, refh, evecs, jstart,
         fconst);
  CHECK(g_den == 1, "write calls den_in_out(1)");
  CHECK(g_mopend == 0, "write no mopend");
  CHECK(std::abs(time - 5.0) < 1e-12, "time>1e7 wrapped to time-1e7");
  {
    std::ifstream f("forsav_test.bin", std::ios::binary);
    double t = 0, r = 0;
    int p = 0, na = 0, no = 0;
    f.read(reinterpret_cast<char*>(&t), sizeof(double));
    f.read(reinterpret_cast<char*>(&p), sizeof(int));
    f.read(reinterpret_cast<char*>(&r), sizeof(double));
    f.read(reinterpret_cast<char*>(&na), sizeof(int));
    f.read(reinterpret_cast<char*>(&no), sizeof(int));
    CHECK(std::abs(t - 5.0) < 1e-12 && p == 2 && std::abs(r + 1.5) < 1e-12,
          "header time=5 ipt=2 refh=-1.5");
    CHECK(na == 1 && no == 4, "header numat=1 norbs=4");
  }

  // Scenario B: read
  std::vector<std::vector<double>> r_dip(4, std::vector<double>(3, -9.0));
  std::vector<double> r_fm(4, -9.0), r_coord(3, -9.0), r_evec(5, -9.0);
  std::vector<double> r_fc(3, -9.0);
  double r_time = -1.0, r_refh = -1.0;
  int r_ipt = 0, r_jstart = -1;
  forsav(r_time, r_dip, r_ipt, r_fm, r_coord, nvar, r_refh, r_evec,
         r_jstart, r_fc);
  CHECK(g_mopend == 0, "read no mopend (valid file)");
  CHECK(g_den == 1, "read does not call den_in_out");
  CHECK(std::abs(r_time - 5.0) < 1e-12 && r_ipt == 2,
        "time/ipt restored 5/2");
  CHECK(std::abs(r_refh + 1.5) < 1e-12, "refh restored -1.5");
  CHECK(std::abs(r_coord[1] - 1.5) < 1e-12 &&
            std::abs(r_coord[2] - 2.5) < 1e-12,
        "coord restored 1.5/2.5");
  CHECK(std::abs(r_fm[3] - 0.7) < 1e-12, "fmatrx(3) restored 0.7");
  CHECK(std::abs(r_dip[2][1] - 0.2) < 1e-12 &&
            std::abs(r_dip[3][2] - 1.3) < 1e-12,
        "deldip restored (2,1)=0.2 (3,2)=1.3");
  CHECK(std::abs(r_evec[4] - 0.6) < 1e-12, "evecs(4) restored 0.6");
  CHECK(r_jstart == 7 && std::abs(r_fc[2] - 3.2) < 1e-12,
        "jstart/fconst restored 7/3.2");

  // Scenario C: missing file -> mopend
  std::remove("forsav_test.bin");
  int g_mopend0 = g_mopend;
  double m_time = 0.0, m_refh = 0.0;
  int m_ipt = 0, m_js = 0;
  std::vector<std::vector<double>> m_dip(4, std::vector<double>(3, 0.0));
  std::vector<double> m_fm(4, 0.0), m_c(3, 0.0), m_ev(5, 0.0), m_fc(3, 0.0);
  forsav(m_time, m_dip, m_ipt, m_fm, m_c, nvar, m_refh, m_ev, m_js, m_fc);
  CHECK(g_mopend == g_mopend0 + 1, "missing file -> mopend once");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
