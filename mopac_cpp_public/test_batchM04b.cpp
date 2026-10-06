// test_batchM04b.cpp  M04 local (Perkins-Stewart LMO localization) batch.
// Validates the Jacobi rotation loop (maximizes sum of (psi^2)^2 over
// atoms), the energy reordering, and the F90 tail behavior: without GRAPH
// the input eigenvectors/energies are restored. resolv/phase_lock/matout/
// molval are unported printing/ambiguity kernels (extern C stubs).
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "local.h"
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

static void setup_local() {
  norbs = 2; numat = 2; nbeta = 0;
  keywrd = " ";
  nfirst.assign({0, 1, 2});
  nlast.assign({0, 1, 2});
  p.assign(4, 0.0);
  pa.assign(4, 0.0);
  pb.assign(4, 0.0);
}

int main() {
  std::printf("M04b batch - local\n");
  double c[4] = {0.8, 0.6, 0.6, -0.8};  // 2x2 col-major-ish (row-major here)
  double eig[2] = {2.0, 1.0};
  double c0[4], eig0[2];
  for (int i = 0; i < 4; ++i) c0[i] = c[i];
  eig0[0] = eig[0]; eig0[1] = eig[1];

  // Test 1: no GRAPH -> input c and eig restored (F90 tail).
  setup_local();
  for (int i = 0; i < 4; ++i) c[i] = c0[i];
  eig[0] = eig0[0]; eig[1] = eig0[1];
  local(c, 2, eig, 0, "c ");
  bool same_c = true, same_eig = true;
  for (int i = 0; i < 4; ++i) if (std::fabs(c[i] - c0[i]) > 1e-12) same_c = false;
  if (std::fabs(eig[0] - eig0[0]) > 1e-12) same_eig = false;
  if (std::fabs(eig[1] - eig0[1]) > 1e-12) same_eig = false;
  CHECK(same_c, "local restores c without GRAPH");
  CHECK(same_eig, "local restores eig without GRAPH");

  // Test 2: GRAPH present -> rotation runs, columns stay normalized.
  setup_local();
  keywrd = " GRAPH";
  for (int i = 0; i < 4; ++i) c[i] = c0[i];
  eig[0] = eig0[0]; eig[1] = eig0[1];
  local(c, 2, eig, 0, "c ");
  bool changed = false;
  for (int i = 0; i < 4; ++i)
    if (std::fabs(c[i] - c0[i]) > 1e-6) changed = true;
  CHECK(changed, "local rotates orbitals with GRAPH");
  double n1 = std::sqrt(c[0] * c[0] + c[1] * c[1]);
  double n2 = std::sqrt(c[2] * c[2] + c[3] * c[3]);
  CHECK(std::fabs(n1 - 1.0) < 1e-9, "local column 1 normalized");
  CHECK(std::fabs(n2 - 1.0) < 1e-9, "local column 2 normalized");
  CHECK(eig[0] <= eig[1] + 1e-12, "local energies reordered ascending");

  // Test 3: already-localized input does no rotation (sum ~ 0 -> no change).
  setup_local();
  for (int i = 0; i < 4; ++i) c[i] = (i % 2 == 0) ? 1.0 : 0.0;  // {1,0,0,1}
  double c1[4];
  for (int i = 0; i < 4; ++i) c1[i] = c[i];
  local(c, 2, eig, 0, "c ");
  bool same2 = true;
  for (int i = 0; i < 4; ++i)
    if (std::fabs(c[i] - c1[i]) > 1e-12) same2 = false;
  CHECK(same2, "local identity on localized input");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
