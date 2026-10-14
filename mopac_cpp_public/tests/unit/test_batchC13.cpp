// test_batchC13.cpp — batch tests for 批C13 (prtdrc).
// T1 first call: init arrays, parref/vref saved; T2 second call runs;
// T3 third call runs the print-decision block; T4 escf<-1e8 early return;
// T5 T-PRI keyword parsing path.
#include "prtdrc.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "drc_C.h"
#include "parameters_C.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace drc_C;
using namespace parameters_C;

static int failures = 0;

static void check(const char* name, bool ok) {
  if (!ok) {
    std::cerr << "FAIL: " << name << "\n";
    ++failures;
  } else {
    std::cout << "ok: " << name << "\n";
  }
}

static void base_setup() {
  numcal = 7;  // fixed distinct calc number so icalcn != numcal triggers init
  moperr = false;
  id = 0;
  keywrd = " ";
  natoms = 1;
  numat = 1;
  nvar = 3;
  nat.assign(2, 1);            // H
  nfirst.assign(2, 1);
  nlast.assign(2, 1);
  na.assign(2, 0);
  nb.assign(2, 0);
  nc.assign(2, 0);
  na_store.assign(2, 0);
  p.assign(4, 0.0);
  std::fill(tore, tore + 2, 0.0);
  tore[1] = 1.0;
  coord.assign(4, std::vector<double>(2, 0.0));
  // itemp_1 (jloop) alias lives in molkst_C
  itemp_1 = 0;
}

static void drive(double escf_v, double ekin_v, double deltt) {
  std::vector<double> xparam(4, 0.0), ref(4, 0.0), velo0(4, 0.0);
  xparam[1] = 0.0; xparam[2] = 0.0; xparam[3] = 0.0;
  ref[1] = 0.0; ref[2] = 0.0; ref[3] = 0.0;
  velo0[1] = 0.1; velo0[2] = 0.0; velo0[3] = 0.0;
  std::vector<std::array<int, 2>> mcoprt;
  double gtot = 0.0, etot = 0.0;
  molkst_C::escf = escf_v;
  prtdrc(deltt, xparam, ref, ekin_v, gtot, etot, velo0, mcoprt, 0, true);
}

int main() {
  base_setup();
  // T1: first call initializes drc_C storage and saves parref/vref
  drive(-100.0, 10.0, 1.0e-15);
  check("T1 vref sized", vref.size() >= 4);
  check("T1 vref saved", std::fabs(vref[1] - 0.1) < 1e-12);
  check("T1 allxyz sized", allxyz.size() >= 10);
  check("T1 parref saved", std::fabs(parref[1] - 0.0) < 1e-12);
  // T2: second call (iloop=2) runs without crash
  drive(-100.0, 10.0, 1.0e-15);
  check("T2 second call ok", !moperr);
  // T3: third call (iloop=3) exercises print-decision block
  drive(-100.0, 10.0, 1.0e-15);
  check("T3 third call ok", !moperr);
  // T4: escf < -1e8 -> early return, storage untouched
  std::size_t sz = vref.size();
  drive(-2.0e9, 0.0, 1.0e-15);
  check("T4 early return keeps vref", vref.size() == sz);
  // T5: T-PRI keyword parsing (stept set), run once more
  base_setup();
  keywrd = " T-PRI=0.05";
  drive(-100.0, 10.0, 1.0e-15);
  check("T5 T-PRI call ok", !moperr);

  if (failures == 0) {
    std::cout << "ALL PASS\n";
    return 0;
  }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
