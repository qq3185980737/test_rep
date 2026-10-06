// test_rotate.cpp — rotate group (rotate.F90 113 lines). Real machine
// assertions + ASan. rotatd_/elenuc_/nddo_to_point_ are deterministic stubs
// (the original Fortran kernels are not part of the 2016 source tree).
#include <cstdio>
#include <cmath>
#include "rotate.h"
#include "parameters_C.h"
#include "molkst_C.h"
#include "rotate_C.h"

// deterministic kernels: enuc=1.0, w block of 3 values 0.5/0.25/0.125,
// kr advances by 4; elenuc fills en(1..6) = -(i*10+k); nddo_to_point scales.
extern "C" void rotatd_(int*, int*, const double*, const double*, double* w,
                        int* kr, double* enuc) {
  *enuc = 1.0;
  w[0] = 0.5; w[1] = 0.25; w[2] = 0.125;
  *kr = 4;
}
extern "C" void elenuc_(int* a1, int* a2, int* a3, int* a4, double* en) {
  // F90 call: elenuc(1, li, li+1, li+lj, en); the args delimit the two
  // atom blocks; the stub fills the whole packed matrix for both atoms.
  int iup = *a4;
  for (int i = 1; i <= iup; ++i)
    for (int k = 1; k <= i; ++k) en[(i*(i-1))/2 + k - 1] = -((double)(i*10 + k));
}
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a,
                               double* enuc, double* rij, int*, int*) {
  for (int k = 0; k < 4; ++k) w[k] *= 2.0;
  *enuc *= 2.0;
}

static int failures = 0;
#define CHECK(cond, msg)                                                     \
  do {                                                                       \
    if (!(cond)) {                                                           \
      std::fprintf(stderr, "FAIL [%s:%d] %s\n", __FILE__, __LINE__, msg);    \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

static void test_small_rij() {
  double xi[3] = {0, 0, 0}, xj[3] = {0.001, 0.001, 0.001}; // rij < 2e-5
  double w[2026] = {0}, e1b[45] = {1,2,3}, e2a[45] = {4,5,6}, enuc = 9.0;
  int kr = 0;
  rotate(1, 2, xi, xj, w, kr, e1b, e2a, enuc);
  CHECK(enuc == 0.0, "small rij: enuc zeroed");
  CHECK(e1b[0] == 0.0 && e2a[0] == 0.0 && w[1] == 0.0, "small rij: arrays zeroed");
}

static void test_normal_path() {
  using namespace parameters_C;
  natorb[1] = 1; natorb[2] = 1;
  double xi[3] = {0, 0, 0}, xj[3] = {1.0, 0, 0};
  double w[2026] = {0}, e1b[45] = {0}, e2a[45] = {0}, enuc = 0.0;
  int kr = 0;
  rotate(1, 2, xi, xj, w, kr, e1b, e2a, enuc);
  CHECK(kr == 4, "rotate: kr from rotatd (4)");
  CHECK(std::fabs(enuc - 1.0) < 1e-12, "rotate: enuc from rotatd");
  CHECK(std::fabs(w[1] - 0.5) < 1e-12 && std::fabs(w[2] - 0.25) < 1e-12 &&
        std::fabs(w[3] - 0.125) < 1e-12, "rotate: w block shifted to 1-based");
  // e1b: li=1 -> ik=1 -> en((1*0)/2+1) = en(1) = -11
  CHECK(std::fabs(e1b[0] + 11.0) < 1e-9, "rotate: e1b(1) = en(1) = -11");
  // e2a: i=2,k=2 -> en((2*1)/2+2)=en(3) = -22
  CHECK(std::fabs(e2a[0] + 22.0) < 1e-9, "rotate: e2a(1) = en(3) = -22");
}

static void test_pm7_to_point() {
  using namespace parameters_C;
  natorb[1] = 1; natorb[2] = 1;
  molkst_C::id = 3; molkst_C::method_pm7 = false;
  double xi[3] = {0, 0, 0}, xj[3] = {1.0, 0, 0};
  double w[2026] = {0}, e1b[45] = {0}, e2a[45] = {0}, enuc = 0.0;
  int kr = 0;
  rotate(1, 2, xi, xj, w, kr, e1b, e2a, enuc);
  CHECK(std::fabs(w[1] - 1.0) < 1e-12 && std::fabs(w[2] - 0.5) < 1e-12,
        "rotate: id==3 routes through nddo_to_point (w doubled)");
  CHECK(std::fabs(enuc - 2.0) < 1e-12, "rotate: id==3 enuc doubled");
  molkst_C::id = 0;
}

int main() {
  test_small_rij();
  test_normal_path();
  test_pm7_to_point();
  if (failures == 0) std::printf("ALL PASS: M05 rotate group\n");
  else std::printf("%d check(s) FAILED\n", failures);
  return failures == 0 ? 0 : 1;
}
