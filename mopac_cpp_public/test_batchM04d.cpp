// test_batchM04d.cpp  M04 mlmo (MOZYME LMO bookkeeping) batch.
// All assertions are hand-derived from mlmo.F90 for the three (ii,jj)
// cases: occupied-only, virtual-only, and both. Padding/expansion of
// cocc/cvir and nf_loc/ne reservation follow the F90 exactly.
#include <cstdio>
#include <vector>

#include "mlmo.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

using namespace molkst_C;
using namespace MOZYME_C;

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

static void setup(int iz[5], int ib[5]) {
  numat = 2; norbs = 4;
  ipad2 = 4; ipad4 = 8;
  for (int i = 0; i <= 4; ++i) iz[i] = ib[i] = 3;
}

int main() {
  std::printf("M04d batch - mlmo\n");
  int iorbs[5] = {0, 2, 2, 0, 0};
  int nce[16] = {0}, ncf[16] = {0}, ncocc[16] = {0}, ncvir[16] = {0};
  int icocc[32] = {0}, icvir[32] = {0};
  std::vector<double> cocc(40, 5.0), cvir(40, 5.0);
  int iz[5], ib[5];

  // Test 1: ii=1, jj=0 (occupied only).
  setup(iz, ib);
  int locc = 1, lvir = 1, nf_loc = 1, ne = 1, nocc = 0, nvir = 0;
  mlmo(locc, lvir, 1, 0, nf_loc, ne, nocc, nvir, iz, ib, nce, ncf, ncocc,
       ncvir, iorbs, icocc, icvir, cocc.data(), cvir.data());
  CHECK(iz[1] == 1, "mlmo occ: iz(1) decremented twice");
  CHECK(ib[1] == 2, "mlmo occ: ib(1) decremented");
  CHECK(nocc == 1 && ncocc[1] == 1, "mlmo occ: ncocc(1)=locc=1");
  CHECK(locc == 1 + 8, "mlmo occ: locc expanded to iocc+k");
  CHECK(nf_loc == 1 + 4, "mlmo occ: nf_loc = nfs+j");
  CHECK(icocc[2] == 1, "mlmo occ: icocc(2)=ii");
  CHECK(ncf[1] == 1, "mlmo occ: ncf(1)=1");
  CHECK(cocc[5] == 0.0, "mlmo occ: cocc padding zeroed");

  // Test 2: ii=0, jj=1 (virtual only).
  setup(iz, ib);
  locc = 1; lvir = 1; nf_loc = 1; ne = 1; nocc = 0; nvir = 0;
  mlmo(locc, lvir, 0, 1, nf_loc, ne, nocc, nvir, iz, ib, nce, ncf, ncocc,
       ncvir, iorbs, icocc, icvir, cocc.data(), cvir.data());
  CHECK(iz[1] == 3, "mlmo vir: iz(1) unchanged (-- then ++)");
  CHECK(ib[1] == 2, "mlmo vir: ib(1) decremented");
  CHECK(nvir == 1 && ncvir[1] == 1, "mlmo vir: ncvir(1)=lvir=1");
  CHECK(ne == 1 + 4, "mlmo vir: ne = nes+j");
  CHECK(nce[1] == 1 && icvir[2] == 1, "mlmo vir: nce/icvir");
  CHECK(lvir == 1 + 8, "mlmo vir: lvir expanded to ivir+k");
  CHECK(cvir[5] == 0.0, "mlmo vir: cvir padding zeroed");

  // Test 3: ii=1, jj=2 (occupied and virtual).
  setup(iz, ib);
  locc = 1; lvir = 1; nf_loc = 1; ne = 1; nocc = 0; nvir = 0;
  mlmo(locc, lvir, 1, 2, nf_loc, ne, nocc, nvir, iz, ib, nce, ncf, ncocc,
       ncvir, iorbs, icocc, icvir, cocc.data(), cvir.data());
  CHECK(iz[1] == 2 && iz[2] == 2, "mlmo both: iz decremented once each");
  CHECK(ib[1] == 2 && ib[2] == 2, "mlmo both: ib decremented");
  CHECK(ncf[1] == 2, "mlmo both: ncf(1)=2");
  CHECK(nce[1] == 2, "mlmo both: nce(1)=2");
  CHECK(icocc[3] == 2, "mlmo both: icocc(3)=jj");
  CHECK(icvir[3] == 2, "mlmo both: icvir(3)=jj");
  CHECK(locc == 1 + 8 && lvir == 1 + 8, "mlmo both: locc/lvir expanded");
  CHECK(nf_loc == 1 + 4 && ne == 1 + 4, "mlmo both: nf_loc/ne reserved");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
