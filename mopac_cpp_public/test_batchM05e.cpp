// test_batchM05e.cpp — molsymy group (molsym/chi/makopr/orient/plato/cartab/symdec)
// Real machine assertions + ASan. External geout/mopend stubbed (not on the
// tested paths; plato's severe-geometry branch is not exercised).
#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include "molsymy.h"
#include "symmetry_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "bldsym.h"
#include "rsp.h"

using namespace symmetry_C;
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

void geout(int) {}
void mopend(const std::string&) {}

static void reset() {
  numat = 0;
  keywrd = " ";
  moperr = 0;
  numcal = 1;
  ielem[0] = 0;
  for (int i = 1; i <= 20; ++i) ielem[i] = 0;
  for (int i = 0; i <= 20; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) elem[j][k][i] = 0.0;
  jelem.assign(21, std::vector<int>(51, 0));
  nat.assign(64, 0);
  nbonds.assign(64, 0);
  ibonds.assign(16, std::vector<int>(64, 0));
  group.assign(100, 0.0);
  jy.assign(7, 0);
  jx.assign(21, std::string());
  nclass = 0;
  nirred = 0;
  name = "";
}

static void test_symdec() {
  reset();
  int ie[21] = {0};
  ie[1] = 1; ie[3] = 1;                      // operations 1 and 3 present
  CHECK(symdec(0, ie) == true, "symdec: n1=0 (C1) always matches");
  CHECK(symdec(5, ie) == true, "symdec: bits {1,3} subset of {1,3}");
  CHECK(symdec(6, ie) == false, "symdec: bit 2 not present -> no match");
  CHECK(symdec(4, ie) == true, "symdec: bit 3 present -> match");
}

static void test_cartab_C1() {
  reset();
  cartab();
  CHECK(name == "C1", "cartab: empty ielem -> C1");
  CHECK(nclass == 1, "cartab: C1 has 1 class");
  CHECK(nirred == 1, "cartab: C1 has 1 irred rep");
}

static void test_cartab_C2() {
  reset();
  ielem[3] = 1;                              // only C2(Z) present
  cartab();
  CHECK(name == "C2", "cartab: ielem[3]=1 -> C2");
  CHECK(nclass == 2, "cartab: C2 has 2 classes");
  CHECK(nirred == 2, "cartab: C2 has 2 irred reps");
  CHECK(jy[2] == 3, "cartab: class 2 is operation 3 (C2z)");
  CHECK(jx[1] == "A", "cartab: irrep 1 name A");
  CHECK(jx[2] == "B", "cartab: irrep 2 name B");
  CHECK(std::fabs(group[(0)*20 + 0] - 1.0) < 1e-12, "cartab: group(1,1)=1");
  CHECK(std::fabs(group[(0)*20 + 1] - 1.0) < 1e-12, "cartab: group(2,1)=1");
  CHECK(std::fabs(group[(1)*20 + 0] - 1.0) < 1e-12, "cartab: group(1,2)=1");
  CHECK(std::fabs(group[(1)*20 + 1] + 1.0) < 1e-12, "cartab: group(2,2)=-1");
}

static void test_chi_C2y() {
  reset();
  numat = 3;
  nat[1] = 7; nat[2] = 6; nat[3] = 7;        // N-C-N linear along Y
  double coords[9] = {0, 1, 0, 0, 0, 0, 0, -1, 0}; // col-major (comp, atom)
  bldsym(2, 2);                               // C2(Y) operation into elem(..,..,2)
  int iqual = 0;
  chi(0.1, coords, 2, iqual);
  CHECK(ielem[2] == 1, "chi: linear along Y is invariant under C2y");
  CHECK(iqual == 3, "chi: all atoms on the C2y axis are fixed points");
  // C2x (perpendicular to the molecular axis) is also a symmetry op of a
  // linear molecule; only the central atom is a fixed point.
  bldsym(1, 1);
  chi(0.1, coords, 1, iqual);
  CHECK(ielem[1] == 1, "chi: linear along Y is invariant under C2x");
  CHECK(iqual == 1, "chi: only the central atom is a fixed point under C2x");
}

static void test_molsym_sphere() {
  reset();
  numat = 1;
  nat[1] = 1;                                 // single atom -> spherical
  double coords[3] = {0, 0, 0};
  double r[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ierror = 0;
  molsym(coords, ierror, r);
  CHECK(ierror == 0, "molsym: single atom no error");
  CHECK(ielem[7] == 1, "molsym: inversion set for sphere");
  CHECK(ielem[8] == 1, "molsym: C3 set for sphere");
  CHECK(ielem[10] == 1, "molsym: C5 set for sphere");
  CHECK(ielem[20] == 1, "molsym: infinite flag set for sphere");
  CHECK(name == "R3", "molsym: single atom -> R3 point group");
}

static void test_molsym_linear() {
  reset();
  numat = 3;
  nat[1] = 7; nat[2] = 6; nat[3] = 7;        // N-C-N linear along Y
  double coords[9] = {0, 1, 0, 0, 0, 0, 0, -1, 0};
  double r[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int ierror = 0;
  molsym(coords, ierror, r);
  CHECK(ierror == 0, "molsym: linear molecule no error");
  CHECK(ielem[20] == 1, "molsym: linear molecule sets infinite flag");
  CHECK(name == "D*h", "molsym: linear N-C-N -> D-infinity-h");
  CHECK(ielem[7] == 1, "molsym: linear N-C-N has inversion (D-inf-h)");
  CHECK(ielem[1] == 1 && ielem[2] == 1 && ielem[3] == 1,
        "molsym: linear N-C-N has C2x/C2y/C2z (D-inf-h)");
}

int main() {
  std::fprintf(stderr, "[t0] symdec\n"); std::fflush(stderr);
  test_symdec();
  std::fprintf(stderr, "[t1] cartab_C1\n"); std::fflush(stderr);
  test_cartab_C1();
  std::fprintf(stderr, "[t2] cartab_C2\n"); std::fflush(stderr);
  test_cartab_C2();
  std::fprintf(stderr, "[t3] chi_C2y\n"); std::fflush(stderr);
  test_chi_C2y();
  std::fprintf(stderr, "[t4] molsym_sphere\n"); std::fflush(stderr);
  test_molsym_sphere();
  std::fprintf(stderr, "[t5] molsym_linear\n"); std::fflush(stderr);
  test_molsym_linear();
  std::fprintf(stderr, "[t6] done\n"); std::fflush(stderr);
  if (failures == 0) {
    std::printf("ALL PASS: M05e molsymy group\n");
    return 0;
  }
  std::fprintf(stderr, "%d check(s) FAILED\n", failures);
  return 1;
}
