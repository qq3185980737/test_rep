// test_hybrid.cpp — hybrid/minloc group. Real machine assertions + ASan.
// External local2/ijbo/rsp are real units (local2.cpp/ijbo.cpp/rsp.cpp).
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include "hybrid.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace molkst_C;

static int failures = 0;
#define CHECK(cond, msg)                                                     \
  do {                                                                       \
    if (!(cond)) {                                                           \
      std::fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);    \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

static void reset() {
  numat = 0;
  keywrd = " ";
  morb = 0;
  norbs = 0;
  iorbs.assign(64, 0);
  nat.assign(64, 0);
  nbonds.assign(64, 0);
  ibonds.assign(16, std::vector<int>(64, 0));
  f.assign(512, 0.0);
  lijbo = true;
  nijbo.assign(64, std::vector<int>(64, 0));
}

// 1) single hydrogen atom: one orbital = the hybrid itself
static void test_h() {
  reset();
  numat = 1;
  morb = 4; norbs = 4;
  nat[1] = 1;
  iorbs[1] = 1;
  nbonds[1] = 0;
  std::vector<double> catom((size_t)morb * norbs, 0.0);
  hybrid(catom.data());
  CHECK(std::fabs(catom[0] - 1.0) < 1e-12, "H: catom(1,1)=1");
  CHECK(catom[1] == 0.0 && catom[2] == 0.0 && catom[3] == 0.0,
        "H: no leakage into other rows/columns");
}

// 2) bare heavy atom Si (9 AOs): 4 sp3 hybrids normalised, 5 d columns identity
static void test_si() {
  reset();
  numat = 1;
  morb = 9; norbs = 9;
  nat[1] = 14;
  iorbs[1] = 9;
  nbonds[1] = 0;
  std::vector<double> catom((size_t)morb * norbs, 0.0);
  hybrid(catom.data());
  for (int col = 1; col <= 4; ++col) {          // sp3 hybrids
    double s = 0.0;
    for (int row = 1; row <= 4; ++row) s += catom[(col-1)*morb + (row-1)] *
                                           catom[(col-1)*morb + (row-1)];
    CHECK(std::fabs(s - 1.0) < 1e-9, "Si: sp3 column normalised");
    for (int row = 5; row <= 9; ++row)
      CHECK(catom[(col-1)*morb + (row-1)] == 0.0, "Si: sp3 column has no d content");
  }
  for (int col = 5; col <= 9; ++col) {          // d columns are the identity set
    for (int row = 1; row <= 9; ++row) {
      double want = (row == col) ? 1.0 : 0.0;
      CHECK(std::fabs(catom[(col-1)*morb + (row-1)] - want) < 1e-12,
            "Si: d column identity");
    }
  }
}

// 3) C-H diatomic: C has 1 bond, ijbo routed through nijbo table,
//    f supplies the overlap block; output columns 1-4 normalised,
//    column 5 (H) is a unit vector.
static void test_ch() {
  reset();
  numat = 2;
  morb = 4; norbs = 5;
  nat[1] = 6; nat[2] = 1;
  iorbs[1] = 4; iorbs[2] = 1;
  nbonds[1] = 1; nbonds[2] = 1;
  ibonds[1][1] = 2;
  ibonds[1][2] = 1;
  // ijbo(2,1) -> nijbo[2][1]; return jl=4 -> f[4+1..4+4] block used (1-based f)
  nijbo[2][1] = 4;
  f[5] = 0.9; f[6] = 0.3; f[7] = 0.3; f[8] = 0.3;   // s-orbital block (m=1)
  std::vector<double> catom((size_t)morb * norbs, 0.0);
  hybrid(catom.data());
  for (int col = 1; col <= 4; ++col) {
    double s = 0.0;
    for (int row = 1; row <= 4; ++row) s += catom[(col-1)*morb + (row-1)] *
                                           catom[(col-1)*morb + (row-1)];
    CHECK(std::fabs(s - 1.0) < 1e-9, "CH: hybrid column normalised");
  }
  // H's column (5): single orbital = 1
  CHECK(std::fabs(catom[(5-1)*morb + 0] - 1.0) < 1e-12, "CH: H column unit");
  for (int row = 2; row <= 4; ++row)
    CHECK(catom[(5-1)*morb + (row-1)] == 0.0, "CH: H column has no leakage");
}

int main() {
  test_h();
  test_si();
  test_ch();
  if (failures == 0) std::printf("ALL PASS: M05 hybrid group\n");
  else std::printf("%d check(s) FAILED\n", failures);
  return failures == 0 ? 0 : 1;
}
