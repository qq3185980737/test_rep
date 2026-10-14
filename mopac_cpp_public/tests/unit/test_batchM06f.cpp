// test_batchM06f.cpp  M06 symtry (dependent internal-coordinate symmetry).
// Exercises haddon across the Cartesian branch (na==0), the angle branch
// (na!=0) and the bond-length multiplier branch (idepfn==19), with exact
// hand-checked values.
#include <cmath>
#include <cstdio>
#include <vector>

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "symtry.h"

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
#define CHK_D(a, b, tol, msg) CHECK(std::fabs((a) - (b)) < (tol), msg)

using namespace common_arrays_C;
using namespace molkst_C;
using namespace symmetry_C;

int main() {
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::printf("M06f batch - symtry\n");

    ndep = 4;
    // na: atoms 1,2 Cartesian (0); atoms 3,4 angle/dihedral type (1).
    na.assign({0, 0, 0, 1, 1});
    geo.assign(5, std::vector<double>(4, 0.0));
    geo[2][1] = 1.2;
    geo[3][1] = 1.5;
    geo[3][3] = 0.5;
    geo[4][1] = 2.0;
    idepfn.assign({0, 1, 4, 8, 19});
    locpar.assign({0, 2, 2, 3, 4});
    locdep.assign({0, 1, 2, 2, 3});
    depmul.assign({0.0, 1.0, 1.0, 1.0, 1.0});
    numcal = 1;

    symtry();

    // i=1: idepfn=1, na[2]==0 -> X=X -> l=1, geo[1][1]=geo[2][1]=1.2.
    CHK_D(geo[1][1], 1.2, 1e-12, "symtry X=X copy");
    // i=2: idepfn=4, na[2]==0 -> X=-X -> geo[1][2]=-geo[2][1]=-1.2.
    CHK_D(geo[1][2], -1.2, 1e-12, "symtry X=-X negate");
    // i=3: idepfn=8, na[3]==1 -> angle branch, w=pi-a(3,3)=pi-0.5, l=3,
    //      stored at geo[3][2].
    CHK_D(geo[3][2], funcon_C::pi - 0.5, 1e-12, "symtry angle pi-a(3)");
    // i=4: idepfn=19, depmul(n+1)=depmul(1)=1>1e-20 -> multiplier branch,
    //      na[4]==1 -> l=1, w=geo[4][1]*fact=2.0, stored at geo[1][3].
    CHK_D(geo[1][3], 2.0, 1e-12, "symtry bond multiplier");

    std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
    return n_fail == 0 ? 0 : 1;
}
