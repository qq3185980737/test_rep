// test_batchC10.cpp — batch tests for 批C10 (hcore_for_MOZYME).
// T1: single H atom, mode=0 fresh build: h diagonal = uspd.  T2: FIELD keyword
// parsing sets efield.  T3: mode=-1 flips h / saves parth & refnuc.  T4: two
// atoms exercise h1elec/rotate/outer paths and the W-packing kr counter.
// T5: debug output path (keywrd "HCORE") runs without crashing.
#include "hcore_for_MOZYME.h"
#include "ijbo.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "chanel_C.h"
#include "funcon_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace cosmo_C;
using namespace chanel_C;
using namespace funcon_C;
using namespace overlaps_C;
using namespace parameters_C;
using namespace common_arrays_C;
using namespace MOZYME_C;

static int cal_ctr = 0;
static int failures = 0;

static void check(const char* name, bool ok) {
  if (!ok) {
    std::cerr << "FAIL: " << name << "\n";
    ++failures;
  } else {
    std::cout << "ok: " << name << "\n";
  }
}

static void set_atom(int i, int z, double x, double y, double zz) {
  nat[i] = z;
  coord[1][i] = x;
  coord[2][i] = y;
  coord[3][i] = zz;
}

static void setup_1atom() {
  numcal = ++cal_ctr;
  moperr = false;
  mozyme = true;
  keywrd = " ";
  numat = 1;
  norbs = 1;
  mpack = 1;
  n2elec = 0;
  id = 0;
  mode = 0;
  useps = false;
  direct = true;
  semidr = true;
  cutofs = 49.0;
  efield[1] = efield[2] = efield[3] = 0.0;
  coord.assign(4, std::vector<double>(20, 0.0));
  nat.assign(20, 0);
  nfirst.assign(20, 0);
  nlast.assign(20, 0);
  l_atom.assign(20, false);   // vecprt_for_MOZYME option-1 path reads l_atom
  iorbs.assign(20, 0);
  jopt.assign(20, 0);
  numred = 1;
  jopt[1] = 1;
  iorbs[1] = 1;
  natorb[1] = 1;
  uspd.assign(2, 0.0);
  uspd[1] = -13.6;
  h.assign(4, 0.0);
  w.assign(2048, 0.0);
  wk.assign(2048, 0.0);
  parth.assign(4, 0.0);
  refnuc = 0.0;
  enuclr = 0.0;
  lijbo = true;
  nijbo.assign(20, std::vector<int>(20, 0));
  nijbo[1][1] = 0;
  iij.assign(20, 0);
  numij.assign(20, 0);
  ijall.assign(20, 0);
  iijj.assign(20, 0);
  tore[1] = 1.0;
  dd[1] = 0.0;
  set_atom(1, 1, 0.0, 0.0, 0.0);
  iw = 6;
}

static void setup_2atom() {
  setup_1atom();
  numat = 2;
  norbs = 2;
  mpack = 3;
  iorbs[2] = 1;
  natorb[2] = 1;
  uspd.assign(3, 0.0);
  uspd[1] = -13.6;
  uspd[2] = -13.6;
  h.assign(4, 0.0);
  parth.assign(4, 0.0);
  nijbo[2][2] = 1;
  nijbo[2][1] = 1;
  nijbo[1][2] = 1;
  tore[1] = 1.0;
  set_atom(2, 1, 1.0, 0.0, 0.0);
}

// T1
static void t_single() {
  setup_1atom();
  mode = 0;
  hcore_for_MOZYME();
  check("T1 h[1]=uspd", std::fabs(h[1] - uspd[1]) < 1e-9);
  check("T1 enuclr=0", std::fabs(enuclr) < 1e-9);
  check("T1 no moperr", !moperr);
}

// T2
static void t_field() {
  setup_1atom();
  mode = 0;
  keywrd = " FIELD=(0.5,0.25,0.0) ";
  hcore_for_MOZYME();
  double expect = 0.5 * a0 / ev;
  check("T2 efield[1]", std::fabs(efield[1] - expect) < 1e-9);
  check("T2 efield[2]", std::fabs(efield[2] - 0.25 * a0 / ev) < 1e-9);
  check("T2 h[1]=uspd", std::fabs(h[1] - uspd[1]) < 1e-9);
}

// T3
static void t_flip() {
  setup_1atom();
  mode = 0;
  hcore_for_MOZYME();
  double h0 = h[1];
  mode = -1;
  hcore_for_MOZYME();
  check("T3 h flipped", std::fabs(h[1] + h0) < 1e-9);
  check("T3 parth saved", std::fabs(parth[1] - h0) < 1e-9);
  check("T3 enuclr flipped", std::fabs(enuclr + 0.0) < 1e-9);
  mode = 1;
  hcore_for_MOZYME();
  check("T3 mode=1 restores h", std::fabs(h[1] - parth[1]) < 1e-9);
  check("T3 refnuc used", std::fabs(enuclr - refnuc) < 1e-9);
}

// T4
static void t_two_atoms() {
  setup_2atom();
  mode = 0;
  hcore_for_MOZYME();
  // h[1] = uspd[1] + e2a*0.5 (from atom 2), h[2] = uspd[2] + e1b*0.5.
  // elenuc stub: en[0]=-1 -> e1b[1]=e2a[1]=-1 -> h += -0.5.
  std::cout << "T4 h[1]=" << h[1] << " uspd=" << uspd[1] << " e2a[0]=" << (h[1]-uspd[1]) << "\n";
  std::cout << "T4 h[2]=" << h[2] << " uspd=" << uspd[2] << " e1b[0]=" << (h[2]-uspd[2]) << "\n";
  check("T4 h[1] = uspd - 1", std::fabs(h[1] - (uspd[1] - 1.0)) < 1e-6);
  check("T4 h[2] = uspd - 1", std::fabs(h[2] - (uspd[2] - 1.0)) < 1e-6);
  check("T4 no moperr", !moperr);
}

// T5
static void t_debug() {
  setup_1atom();
  mode = 0;
  keywrd = " HCORE ";
  l_atom[1] = true;  // vecprt option-1 prints selected atoms; all-false would yield numb=0 (infinite loop in print)
  hcore_for_MOZYME();
  check("T5 debug no moperr", !moperr);
}

int main() {
  std::cout << "start\n"; std::cout.flush();
  t_single();
  std::cout << "after T1\n"; std::cout.flush();
  t_field();
  std::cout << "after T2\n"; std::cout.flush();
  t_flip();
  std::cout << "after T3\n"; std::cout.flush();
  t_two_atoms();
  std::cout << "after T4\n"; std::cout.flush();
  t_debug();
  std::cout << "after T5\n"; std::cout.flush();
  if (failures == 0) {
    std::cout << "ALL PASS\n";
    return 0;
  }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
