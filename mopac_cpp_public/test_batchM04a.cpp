// test_batchM04a.cpp  M04 isitsc (MOZYME SCF convergence test) batch.
// Hand-derived branch checks: convergence via ovmax/energy_diff, rapid exit
// when energies are consistently lower/higher than the previous minimum,
// and the niter>itrmax fallback.
#include <cmath>
#include <cstdio>

#include "isitsc.h"
#include "molkst_C.h"

using namespace molkst_C;
extern double ovmax;   // plain global (MOZYME_C.cpp)
extern double energy_diff;

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
  std::printf("M04a batch - isitsc\n");
  double selcon = 0.01;

  // Test 1: two successive "close" iterations -> converged (iscf=1).
  {
    int iemin = 0, iemax = 0;
    bool okscf = false;
    ovmax = 1e-6;      // < fmo_test = 0.05
    energy_diff = 1e-6;  // < energy_test = 0.01
    isitsc(0.5, selcon, 0.4, iemin, iemax, okscf, 1, 50);
    CHECK(!okscf, "isitsc first close pass not converged (scf1 priming)");
    isitsc(0.5, selcon, 0.4, iemin, iemax, okscf, 2, 50);
    CHECK(okscf, "isitsc second close pass converged");
    CHECK(iscf == 1, "isitsc iscf=1 for scf1 path");
  }

  // Test 2: niter > itrmax forces convergence even with a large gap.
  {
    int iemin = 0, iemax = 0;
    bool okscf = false;
    ovmax = 1.0;
    energy_diff = 1.0;
    isitsc(0.5, selcon, 0.4, iemin, iemax, okscf, 60, 50);
    CHECK(okscf, "isitsc niter>itrmax fallback");
  }

  // Test 3: consistently lower energies -> rapid exit iscf=1.
  {
    int iemin = 0, iemax = 0;
    bool okscf = false;
    ovmax = 1.0;  // keep ovmax far above fmo_test so scf1=false
    energy_diff = 1.0;
    double emin = 10.0;
    // Successive deltas 0.01/0.01/0.009 must stay under
    // 0.1*(emin-escf_current)=0.0129 at the final call.
    double escf = 9.90;
    isitsc(escf, selcon, emin, iemin, iemax, okscf, 1, 50);
    CHECK(!okscf, "isitsc lower-energy run 1 not yet converged");
    escf = 9.89;
    isitsc(escf, selcon, emin, iemin, iemax, okscf, 2, 50);
    CHECK(!okscf, "isitsc lower-energy run 2 not yet converged");
    escf = 9.88;
    isitsc(escf, selcon, emin, iemin, iemax, okscf, 3, 50);
    CHECK(!okscf, "isitsc lower-energy run 3 not yet converged");
    escf = 9.871;
    isitsc(escf, selcon, emin, iemin, iemax, okscf, 4, 50);
    CHECK(okscf, "isitsc lower-energy rapid exit");
    CHECK(iscf == 1, "isitsc rapid-exit iscf=1");
  }

  // Test 4: a large energy jump breaks the rapid-exit chain.
  {
    int iemin = 0, iemax = 0;
    bool okscf = false;
    ovmax = 1.0;
    energy_diff = 1.0;
    double emin = 1.0;
    isitsc(0.95, selcon, emin, iemin, iemax, okscf, 1, 50);
    isitsc(0.94, selcon, emin, iemin, iemax, okscf, 2, 50);
    isitsc(0.93, selcon, emin, iemin, iemax, okscf, 3, 50);
    isitsc(0.5, selcon, emin, iemin, iemax, okscf, 4, 50);
    CHECK(!okscf, "isitsc jump breaks rapid exit");
  }

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
