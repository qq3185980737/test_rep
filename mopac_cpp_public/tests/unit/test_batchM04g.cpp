// test_batchM04g.cpp  M04 lyse (break inter-protein bonds) batch.
// Scenes: H atom keeps only its shortest bond (>0.95A), an O-S bond is
// dropped and the O keeps its original bond count (F90 cycle semantics),
// a C-N bond survives intact, and an O attached to a sulfate (nbonds==4)
// is left un-lysed.
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "lyse.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

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

extern "C" double distance_(int* a, int* b) {
  // Atom 2 at 1.5 qualifies; atom 3 at 0.9 is below the 0.95 cutoff.
  return (*a == 2) ? 1.5 : 0.9;
}

static void setup_arrays() {
  nat.assign(8, 0);
  nbonds.assign(8, 0);
  ibonds.assign(8, std::vector<int>(8, 0));
}

int main() {
  std::printf("M04g batch - lyse\n");

  // Test 1: H(1) bonded to C(2) and O(3); keep only the shortest (atom 2).
  {
    setup_arrays();
    numat = 3;
    nat[1] = 1; nat[2] = 6; nat[3] = 8;
    nbonds[1] = 2; ibonds[1][1] = 2; ibonds[2][1] = 3;
    nbonds[2] = 1; ibonds[1][2] = 1;
    nbonds[3] = 1; ibonds[1][3] = 1;
    lyse();
    CHECK(nbonds[1] == 1 && ibonds[1][1] == 2, "lyse H keeps shortest bond");
    CHECK(nbonds[2] == 1 && ibonds[1][2] == 1, "lyse C-H bond rebuilt");
    CHECK(nbonds[3] == 0, "lyse O bond removed");
  }

  // Test 2: O(1)-S(2): S bond dropped, no N/C neighbor -> nbonds(1) keeps
  // its original value (F90 cycle semantics, no nbonds(i)=l write).
  {
    setup_arrays();
    numat = 2;
    nat[1] = 8; nat[2] = 16;
    nbonds[1] = 1; ibonds[1][1] = 2;
    nbonds[2] = 1; ibonds[1][2] = 1;
    lyse();
    CHECK(nbonds[1] == 1, "lyse O-S: nbonds(1) untouched (cycle)");
  }

  // Test 3: C(1)-N(2) protein bond survives (nat=6 copy path).
  {
    setup_arrays();
    numat = 2;
    nat[1] = 6; nat[2] = 7;
    nbonds[1] = 1; ibonds[1][1] = 2;
    nbonds[2] = 1; ibonds[1][2] = 1;
    lyse();
    CHECK(nbonds[1] == 1 && ibonds[1][1] == 2, "lyse C-N intact");
  }

  // Test 4: O attached to sulfate (S with 4 bonds) is left un-lysed.
  {
    setup_arrays();
    numat = 3;
    nat[1] = 8; nat[2] = 16; nat[3] = 6;
    nbonds[1] = 2; ibonds[1][1] = 2; ibonds[2][1] = 3;
    nbonds[2] = 4; ibonds[1][2] = 1; ibonds[2][2] = 3;
    nbonds[3] = 1; ibonds[1][3] = 1;
    lyse();
    CHECK(nbonds[1] == 2 && ibonds[1][1] == 2, "lyse sulfate O un-lysed");
  }

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
