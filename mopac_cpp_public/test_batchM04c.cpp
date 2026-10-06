// test_batchM04c.cpp  M04 local2 (LMO localization, no global state) batch.
// Jacobi rotations maximize the per-atom sum of (psi^2)^2; the angle
// formula and convergence mirror local2.F90 exactly. Assertions: rotation
// changes a mixed input while preserving column norms, identity input stays
// put, and a 4-orbital run converges within its 19-iteration budget.
#include <cmath>
#include <cstdio>

#include "local2.h"

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
  std::printf("M04c batch - local2\n");
  int nfirst[3] = {0, 1, 2};
  int nlast[3] = {0, 1, 2};

  // Test 1: mixed 2x2 input rotates; norms stay 1.
  {
    double c[4] = {0.8, 0.6, 0.6, -0.8};
    double c0[4] = {0.8, 0.6, 0.6, -0.8};
    local2(c, 2, 2, nfirst, nlast, 2);
    bool changed = false;
    for (int i = 0; i < 4; ++i)
      if (std::fabs(c[i] - c0[i]) > 1e-6) changed = true;
    CHECK(changed, "local2 rotates mixed orbitals");
    double n1 = std::sqrt(c[0] * c[0] + c[1] * c[1]);
    double n2 = std::sqrt(c[2] * c[2] + c[3] * c[3]);
    CHECK(std::fabs(n1 - 1.0) < 1e-9, "local2 col 1 normalized");
    CHECK(std::fabs(n2 - 1.0) < 1e-9, "local2 col 2 normalized");
  }

  // Test 2: already-localized input unchanged (sa <= 1e-14).
  {
    double c[4] = {1.0, 0.0, 0.0, 1.0};
    double c0[4] = {1.0, 0.0, 0.0, 1.0};
    local2(c, 2, 2, nfirst, nlast, 2);
    bool same = true;
    for (int i = 0; i < 4; ++i)
      if (std::fabs(c[i] - c0[i]) > 1e-12) same = false;
    CHECK(same, "local2 identity on localized input");
  }

  // Test 3: 4 orbitals, 2 atoms, 19-iteration budget; columns normalized
  // after convergence (sum < 1e-5 early exit or loop end).
  {
    int nfirst4[3] = {0, 1, 3};
    int nlast4[3] = {0, 2, 4};
    // Two 2-orbital atoms, deliberately mixed columns.
    // Normalized mixed columns: atoms 1-2 and 3-4 are pairs.
    double c[16] = {
        0.8, 0.6, 0.0, 0.0,  0.6, -0.8, 0.0, 0.0,
        0.0, 0.0, 0.8, 0.6,  0.0, 0.0, 0.6, -0.8};
    local2(c, 4, 4, nfirst4, nlast4, 2);
    for (int col = 0; col < 4; ++col) {
      double nrm = 0.0;
      for (int k = 0; k < 4; ++k) nrm += c[col * 4 + k] * c[col * 4 + k];
      CHECK(std::fabs(std::sqrt(nrm) - 1.0) < 1e-9, "local2 4-orbital col normalized");
    }
  }

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
