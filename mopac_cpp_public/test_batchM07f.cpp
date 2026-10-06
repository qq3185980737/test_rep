// test_batchM07f.cpp  M07 drcout (DRC output) batch.
// Non-DRC path: escf3={1,.5,.1}, ekin3={2,.2,.01}, etot3={3,.7,.11},
// xtot3={0,1,0}, fract=0.5 -> rc_escf=1.275, ekin=2.1025, etot=3.3775,
// errr=0 (clamped), rxn_coord=0.5. 2nd call same fract is skipped by the
// rxn_coord duplicate check (jloop rolls back). 3rd call fract=0.7 prints.
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "drcout.h"
#include "elemts_C.h"
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

// Test stubs for external hooks.
static int g_traj = 0;
extern "C" void to_screen_(const char*) {}
void local(double*, int, double*, int, const char*) {}
void mullik() {}
void l_control(const char*, int, int) {}
void write_trajectory(double*, int, double*, double, double, double,
                      double) {
  ++g_traj;
}
double reada(const std::string&, int) { return 0.0; }

int main() {
  std::printf("M07f batch - drcout\n");
  using namespace maps_C;
  using namespace molkst_C;
  using namespace common_arrays_C;

  numcal = 1;
  keywrd = " ";
  natoms = 2;
  numat = 2;
  nat.resize(3, 0);
  nat[1] = 6; nat[2] = 6;
  labels.resize(3, 0);
  labels[1] = 6; labels[2] = 6;
  na.resize(3, 0); nb.resize(3, 0); nc.resize(3, 0);
  loc.resize(3);
  loc[1].assign(2, 0); loc[2].assign(2, 0);
  c.resize(2); c[1].assign(2, 0.0);
  eigs.assign(2, 0.0);

  std::vector<std::vector<double>> xyz3(4, std::vector<double>(7, 0.0));
  std::vector<std::vector<double>> geo3(4, std::vector<double>(7, 0.0));
  std::vector<std::vector<double>> vel3(4, std::vector<double>(7, 0.0));
  geo3[1][1] = 1.0; geo3[1][4] = 2.0; geo3[2][1] = 0.5; geo3[3][1] = 0.1;
  std::vector<double> escf3(4, 0.0), ekin3(4, 0.0), etot3(4, 0.0),
      xtot3(4, 0.0);
  escf3[1] = 1.0; escf3[2] = 0.5; escf3[3] = 0.1;
  ekin3[1] = 2.0; ekin3[2] = 0.2; ekin3[3] = 0.01;
  etot3[1] = 3.0; etot3[2] = 0.7; etot3[3] = 0.11;
  xtot3[1] = 0.0; xtot3[2] = 1.0; xtot3[3] = 0.0;
  std::vector<double> charge(3, 0.0);
  int jloop = 0;

  drcout(xyz3, geo3, vel3, 6, 0.0, escf3, ekin3, etot3, xtot3, 2, charge,
         0.5, " ", " ", 0, jloop);
  CHK_D(rc_escf, 1.275, 1e-9, "drcout rc_escf=1.275");
  CHK_D(ekin, 2.1025, 1e-9, "drcout ekin=2.1025");
  CHK_D(rxn_coord, 0.5, 1e-9, "drcout rxn_coord=0.5");
  CHECK(jloop == 1, "drcout jloop=1 after first call");
  CHECK(line.find("1.27500") != std::string::npos, "drcout line has escf");
  CHECK(line.find("2.10250") != std::string::npos, "drcout line has ekin");
  CHECK(line.find("3.37750") != std::string::npos, "drcout line has total");
  CHECK(g_traj == 1, "drcout trajectory written once");

  // Same fract -> duplicate rxn_coord: roll back and return.
  drcout(xyz3, geo3, vel3, 6, 0.0, escf3, ekin3, etot3, xtot3, 2, charge,
         0.5, " ", " ", 0, jloop);
  CHECK(jloop == 1, "drcout duplicate point rolls jloop back to 1");
  CHECK(g_traj == 1, "drcout no trajectory on duplicate");

  // New fract -> prints again.
  drcout(xyz3, geo3, vel3, 6, 0.0, escf3, ekin3, etot3, xtot3, 2, charge,
         0.7, " ", " ", 0, jloop);
  CHK_D(rxn_coord, 0.7, 1e-9, "drcout rxn_coord=0.7");
  CHECK(jloop == 2, "drcout jloop=2 after 3rd call");
  CHECK(line.find("0.7000") != std::string::npos, "drcout line has 0.7 coord");
  CHECK(g_traj == 2, "drcout trajectory written twice");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
