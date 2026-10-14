// test_batchM05a.cpp  M05a: xyzint/bangle/dihed/dang/dist2/dot1.
// Cartesian<->internal coordinate primitives, all id=0 (plus periodic
// equivalence check for bangle with zero lattice vectors).
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "bangle.h"
#include "dang.h"
#include "dihed.h"
#include "dist2.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "xyzint.h"

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

static const double PI = 3.14159265358979323846;

int main() {
  std::printf("M05a batch - xyzint/bangle/dihed/dang/dist2/dot1\n");
  std::fflush(stdout);
  molkst_C::id = 0;
  molkst_C::l11 = 0;
  molkst_C::l21 = 0;
  molkst_C::l31 = 0;
  common_arrays_C::tvec.assign(4, std::vector<double>(4, 0.0));

  // xyz(3, n): 1-based [d][atom].
  std::vector<std::vector<double>> xyz(4, std::vector<double>(5, 0.0));
  // Atom 1 (1,0,0), atom 2 (0,0,0), atom 3 (0,1,0), atom 4 (0,1,0) planar.
  xyz[1][1] = 1.0; xyz[2][1] = 0.0; xyz[3][1] = 0.0;
  xyz[1][2] = 0.0; xyz[2][2] = 0.0; xyz[3][2] = 0.0;
  xyz[1][3] = 0.0; xyz[2][3] = 1.0; xyz[3][3] = 0.0;
  xyz[1][4] = 0.0; xyz[2][4] = 1.0; xyz[3][4] = 0.0;

  // --- bangle: right angle at atom 2 ---
  double ang = 0.0;
  bangle(xyz, 1, 2, 3, ang);
  std::printf("  bangle(1,2,3) = %.10f\n", ang);
  CHECK(std::abs(ang - PI / 2.0) < 1e-10, "bangle right angle = pi/2");

  // --- bangle periodic branch with zero lattice vectors == id=0 result ---
  molkst_C::id = 1;
  double ang2 = 0.0;
  bangle(xyz, 1, 2, 3, ang2);
  molkst_C::id = 0;
  CHECK(std::abs(ang2 - ang) < 1e-12,
        "bangle periodic branch (zero tvec) matches id=0");

  // --- dang: signed angle ---
  double a1 = 1.0, a2 = 0.0, b1 = 0.0, b2 = 1.0, rcos = 0.0;
  dang(a1, a2, b1, b2, rcos);
  std::printf("  dang = %.10f\n", rcos);
  // F90: sinth>0 -> rcos = -(2pi - acos) = -4.712389
  CHECK(std::abs(rcos - (-4.7123889803846898577)) < 1e-10,
        "dang (1,0)-(0,1) = -3pi/2");

  // --- dihed: planar quadrilateral -> 0 ---
  double dih = 99.0;
  dihed(xyz, 1, 2, 3, 4, dih);
  std::printf("  dihed planar = %.10f\n", dih);
  CHECK(std::abs(dih) < 1e-10, "dihed planar = 0");

  // --- dihed: 45-degree twist ---
  // j=(0,0,0), k=(1,0,0) axis (x). i=(1,1,0): y-proj at 0 deg.
  // l=(1,1,1): k->l=(0,1,1): yz-proj at 45 deg -> dihedral = pi/4.
  std::vector<std::vector<double>> xyz2(4, std::vector<double>(5, 0.0));
  xyz2[1][1] = 1.0; xyz2[2][1] = 1.0; xyz2[3][1] = 0.0;  // i
  xyz2[1][2] = 0.0; xyz2[2][2] = 0.0; xyz2[3][2] = 0.0;  // j
  xyz2[1][3] = 1.0; xyz2[2][3] = 0.0; xyz2[3][3] = 0.0;  // k
  xyz2[1][4] = 1.0; xyz2[2][4] = 1.0; xyz2[3][4] = 1.0;  // l
  double dih2 = 99.0;
  dihed(xyz2, 1, 2, 3, 4, dih2);
  std::printf("  dihed 45deg = %.10f\n", dih2);
  CHECK(std::abs(dih2 - PI / 4.0) < 1e-9, "dihed 45-degree twist = pi/4");

  // --- dist2 / dot1 ---
  std::array<double, 3> p = {1.0, 2.0, 3.0}, q = {1.0, 0.0, 3.0};
  CHECK(std::abs(dist2(p, q) - 4.0) < 1e-15, "dist2 = 4");
  CHECK(std::abs(dot1(p, q) - 10.0) < 1e-15, "dot1 = 10");

  // --- xyzint: 3-atom bent molecule ---
  // O(0,0,0), H1(1,0,0), H2(0,1,0): bond lengths 1, angle 90 deg.
  // geo column-major: geo(d,i) at (i-1)*3+(d-1); 1-based direct indexing.
  std::vector<double> xyzf(9, 0.0);
  xyzf[0] = 0.0; xyzf[1] = 0.0; xyzf[2] = 0.0;  // atom1 O
  xyzf[3] = 1.0; xyzf[4] = 0.0; xyzf[5] = 0.0;  // atom2 H1
  xyzf[6] = 0.0; xyzf[7] = 1.0; xyzf[8] = 0.0;  // atom3 H2
  std::vector<int> na(4, 0), nb(4, 0), nc(4, 0);
  std::vector<double> geo(9, 99.0);
  xyzint(xyzf.data(), 3, na.data(), nb.data(), nc.data(), 57.29578, geo.data());
  // na/nb/nc keep Fortran 1-based direct indexing: na(2) lives in na[2].
  std::printf("  xyzint geo: r(O-H1)=%.6f r(O-H2)=%.6f a=%.6f na2=%d\n",
              geo[3], geo[6], geo[7], na[2]);
  CHECK(std::abs(geo[3] - 1.0) < 1e-10, "xyzint r(O-H1) = 1");
  CHECK(std::abs(geo[6] - 1.0) < 1e-10, "xyzint r(O-H2) = 1");
  CHECK(std::abs(geo[7] - 90.0) < 1e-3,
        "xyzint angle H1-O-H2 = 90 deg (+/- degree-const rounding)");
  CHECK(na[1] == 0, "xyzint na(1)=0 after cleanup");
  CHECK(na[2] == 1, "xyzint na(2)=1 after cleanup");

  std::printf("RESULT: %d PASS, %d FAIL\n", n_pass, n_fail);
  return n_fail == 0 ? 0 : 1;
}
