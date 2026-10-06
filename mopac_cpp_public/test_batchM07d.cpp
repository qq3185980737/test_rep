// test_batchM07d.cpp  M07 picopt (SCF atom selection) batch.
// Scenarios: loop=-1 selects all real atoms; loop=0 selects atoms marked
// in loc plus symmetry-dependent (locdep); on repeated call with same
// numcal, neighbors of optimized atoms are marked (iopt=1).
#include <cstdio>
#include <vector>

#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "picopt.h"
#include "symmetry_C.h"

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

int main() {
  std::printf("M07d batch - picopt\n");
  using common_arrays_C::loc;
  using common_arrays_C::nbonds;
  using common_arrays_C::ibonds;
  using common_arrays_C::labels;
  using molkst_C::natoms;
  using molkst_C::numat;
  using molkst_C::nvar;
  using molkst_C::ndep;
  using molkst_C::numcal;
  using symmetry_C::locdep;

  extern std::vector<int> jopt;
  extern int numred;
  jopt.resize(10);
  natoms = 3; numat = 3; nvar = 1; ndep = 0; numcal = 1;
  labels.resize(4); labels[1] = 6; labels[2] = 6; labels[3] = 99;
  loc.resize(3); loc[1].assign(2, 0); loc[1][1] = 1;
  locdep.resize(2, 0);
  nbonds.resize(4, 0); ibonds.resize(4);
  for (int i = 1; i <= 3; ++i) ibonds[i].assign(3, 0);
  nbonds[1] = 1; ibonds[1][1] = 2;

  // First call: imol(0) != numcal(1), no neighbor propagation.
  numred = 0;
  picopt(0);
  CHECK(numred == 1, "picopt first call selects 1 atom");
  CHECK(jopt[1] == 1, "picopt jopt(1)=1");

  // Second call with same numcal: neighbors of iopt==2 marked.
  picopt(0);
  CHECK(numred == 2, "picopt 2nd call propagates neighbor");
  CHECK(jopt[1] == 1 && jopt[2] == 2,
        "picopt jopt={1,2}");

  // Symmetry-dependent atom 2 (no loc atoms this time).
  nvar = 0; ndep = 1; locdep[1] = 2; nbonds[2] = 0;
  picopt(0);
  CHECK(numred == 1, "picopt locdep selects atom 2");
  CHECK(jopt[1] == 2, "picopt jopt(1)=2");

  // loop=-1: all real atoms.
  picopt(-1);
  CHECK(numred == 3, "picopt loop=-1 all atoms");
  CHECK(jopt[1] == 1 && jopt[2] == 2 &&
            jopt[3] == 3,
        "picopt jopt={1,2,3}");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
