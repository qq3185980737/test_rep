// test_batchM07i.cpp  M07 ef_C/drc_C/derivs_C module globals batch.
// Verifies the three Fortran module mappings (ef_C.F90 1035B: ef_C,
// drc_C, derivs_C) expose every declared symbol and default-init state,
// and that cross-module writes are visible (linkage intact).
#include <cmath>
#include <cstdio>
#include <vector>

#include "derivs_C.h"
#include "drc_C.h"
#include "ef_C.h"

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
  std::printf("M07i batch - ef_C/drc_C/derivs_C modules\n");
  // ef_C scalars (Fortran: ef_mode,nstep,negreq,iprnt,iloop int;
  // ddx,rmin,rmax,omin,xlamd,xlamd0,skal,x0,x1,x2 double).
  CHECK(ef_C::nstep == 0 && ef_C::negreq == 0 && ef_C::iprnt == 0,
        "ef_C ints default 0");
  CHECK(ef_C::ef_mode == 0 && ef_C::iloop == 0, "ef_C ef_mode/iloop 0");
  CHECK(std::abs(ef_C::ddx - 1.0) < 1e-15, "ef_C ddx default 1.0");
  CHECK(ef_C::xlamd == 0.0 && ef_C::xlamd0 == 0.0 && ef_C::skal == 1.0,
        "ef_C xlamd/xlamd0 0 skal 1");
  CHECK(ef_C::rmin == 0.0 && ef_C::rmax == 0.0 && ef_C::omin == 0.0,
        "ef_C rmin/rmax/omin 0");
  CHECK(ef_C::x0 == 0.0 && ef_C::x1 == 0.0 && ef_C::x2 == 0.0,
        "ef_C x0/x1/x2 0");
  // ef_C arrays default empty; alparm resizable 1-based
  CHECK(ef_C::hess.empty() && ef_C::bmat.empty() && ef_C::u.empty(),
        "ef_C hess/bmat/u empty at start");
  CHECK(ef_C::oldhss.empty() && ef_C::oldu.empty() && ef_C::pmat.empty(),
        "ef_C oldhss/oldu/pmat empty");
  CHECK(ef_C::uc.empty() && ef_C::hessc.empty() && ef_C::oldf.empty() &&
            ef_C::d.empty() && ef_C::vmode.empty(),
        "ef_C uc/hessc/oldf/d/vmode empty");
  ef_C::alparm.assign(4, std::vector<double>(5, 7.5));
  CHECK(std::abs(ef_C::alparm[3][4] - 7.5) < 1e-12,
        "ef_C alparm 1-based write/read");
  // drc_C
  CHECK(drc_C::vref.empty() && drc_C::vref0.empty(),
        "drc_C vref/vref0 empty");
  CHECK(drc_C::allxyz.empty() && drc_C::allvel.empty(),
        "drc_C allxyz/allvel empty");
  CHECK(drc_C::xyz3.empty() && drc_C::vel3.empty(),
        "drc_C xyz3/vel3 empty");
  CHECK(drc_C::allgeo.empty() && drc_C::geo3.empty(),
        "drc_C allgeo/geo3 empty");
  CHECK(drc_C::parref.empty() && drc_C::time == 0.0,
        "drc_C parref empty time 0");
  drc_C::time = 42.0;
  drc_C::vref.assign(3, 1.0);
  CHECK(std::abs(drc_C::time - 42.0) < 1e-12 && drc_C::vref[2] == 1.0,
        "drc_C time/vref write visible");
  // derivs_C
  CHECK(derivs_C::wmat.empty() && derivs_C::hmat.empty() &&
            derivs_C::fmat.empty(),
        "derivs_C wmat/hmat/fmat empty");
  CHECK(derivs_C::b.empty() && derivs_C::ab.empty() &&
            derivs_C::fb_ci.empty(),
        "derivs_C b/ab/fb_ci empty");
  CHECK(derivs_C::aidref.empty() && derivs_C::work2.empty(),
        "derivs_C aidref/work2 empty");
  derivs_C::b.assign(4, 2.5);
  CHECK(std::abs(derivs_C::b[3] - 2.5) < 1e-12,
        "derivs_C b 1-based write/read");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
