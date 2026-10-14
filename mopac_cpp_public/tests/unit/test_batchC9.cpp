// test_batchC9.cpp — batch tests for 批C9 (new_esp).
// T1: rys quadrature moments.  T2: vint one-dimensional integrals.
// T3: sqrdc/sqrsl least-squares.  T4: density_for_GPU.  T5: new_esp end-to-end
// (1-atom C, CUBE keyword, writes test_c9.grd).
#include "new_esp.h"
#include "esp_utilities.h"
#include "sqrdc_sqrsl.h"
#include "esp_support.h"
#include "mult.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
#include "overlaps_C.h"
#include "funcon_C.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace overlaps_C;
using namespace parameters_C;

static int cal_ctr = 0;
static int failures = 0;

static void setup() {
  numcal = ++cal_ctr;
  moperr = false;
  mozyme = false;
  keywrd = " ";
  jobnam = "test_c9";
  numat = 0;
  norbs = 0;
  nclose = 0;
  nopen = 0;
  fract = 1.0;
  coord.assign(4, std::vector<double>(20, 0.0));
  nat.assign(20, 0);
  nfirst.assign(20, 0);
  nlast.assign(20, 0);
  c.assign(4, std::vector<double>(20, 0.0));
  p.assign(2, 0.0);
}

static void check(const char* name, bool ok) {
  if (!ok) {
    std::cerr << "FAIL: " << name << "\n";
    ++failures;
  } else {
    std::cout << "ok: " << name << "\n";
  }
}

// T1
static void t_rys() {
  double u[4] = {0, 0, 0, 0}, w[4] = {0, 0, 0, 0};
  rys(1.0, 1, u, w);
  double I0 = 0.5 * std::sqrt(3.141592653589793) * std::erf(1.0);
  check("T1 w1=I0", std::fabs(w[1] - I0) < 1e-5);
  check("T1 u1 in (0,1)", u[1] > 0.0 && u[1] < 1.0);
  // 2-root: sum of weights = I0, weighted u^2 = I1.
  rys(2.0, 2, u, w);
  double I0b = 0.5 * std::sqrt(3.141592653589793) * std::erf(std::sqrt(2.0)) / std::sqrt(2.0);
  double I1b = (I0b - std::exp(-2.0)) / 4.0;
  check("T1 w1+w2=I0", std::fabs(w[1] + w[2] - I0b) < 1e-4);
  check("T1 u^2 weighted", std::fabs(w[1] * u[1] * u[1] + w[2] * u[2] * u[2] - I1b) < 1e-4);
}

// T2
static void t_vint() {
  double xi, yi, zi;
  vint(xi, yi, zi, 0, 0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 0.7);
  check("T2 k=l=0 -> 1", std::fabs(xi - 1.0) < 1e-9 && std::fabs(yi - 1.0) < 1e-9 &&
                            std::fabs(zi - 1.0) < 1e-9);
  // k=1, l=0: mean of linear factor over the Hermite density = x0.
  vint(xi, yi, zi, 1, 0, 1.5, 2.5, 3.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.8);
  check("T2 k=1,l=0 -> x0", std::fabs(xi - 1.5) < 1e-6);
}

// T3
static void t_qrls() {
  // A = [[1,0],[0,1],[1,1]], y = [1,2,3]; least squares x = [1,2].
  int ldx = 3, n = 3, p = 2;
  std::vector<float> x(ldx * p + 1, 0.0f);
  // column-major 1-based: x[(j-1)*ldx + (i-1)] = A(i,j)
  // 1-based array: element (i,j) at x[(j-1)*ldx+i]
  x[(0) * ldx + 1] = 1.0f; x[(0) * ldx + 2] = 0.0f; x[(0) * ldx + 3] = 1.0f;  // col1
  x[(1) * ldx + 1] = 0.0f; x[(1) * ldx + 2] = 1.0f; x[(1) * ldx + 3] = 1.0f;  // col2
  std::vector<float> qraux(p + 1, 0.0f), work(p + 1, 0.0f);
  std::vector<int> jpvt(p + 1, 0);
  sqrdc(x, ldx, n, p, qraux, jpvt, work, 0);
  std::vector<float> y = {0, 1, 2, 3};
  std::vector<float> qy(n + 1, 0.0f), qty(n + 1, 0.0f), b(p + 1, 0.0f);
  std::vector<float> rsd(n + 1, 0.0f), xb(p + 1, 0.0f);
  int info = 0;
  sqrsl(x, ldx, n, p, qraux, y, qy, qty, b, rsd, xb, 111, info);
  check("T3 xb[1]~1", std::fabs(xb[1] - 1.0f) < 1e-3f);
  check("T3 xb[2]~2", std::fabs(xb[2] - 2.0f) < 1e-3f);
  // residual = y - A x = [1-1, 2-2, 3-3] = [0,0,0]
  check("T3 rsd small", std::fabs(rsd[1]) < 1e-3f && std::fabs(rsd[2]) < 1e-3f &&
                           std::fabs(rsd[3]) < 1e-3f);
}

// T4
static void t_density() {
  // C = I (2 AO, 1 MO closed): p(1,1)=2, p(2,1)=0, p(2,2)=2.
  double vecs[4] = {1.0, 0.0, 0.0, 1.0};  // row-major C(i,k)
  std::vector<double> p(4, 0.0);
  density_for_GPU(vecs, 1.0, 1, 1, 2.0, 3, 2, 2, p, 5);
  check("T4 p11=2", std::fabs(p[1] - 2.0) < 1e-12);
  check("T4 p21=0", std::fabs(p[2]) < 1e-12);
  check("T4 p22=0", std::fabs(p[3]) < 1e-12);
}

// T5
static void t_new_esp() {
  numat = 1;
  norbs = 1;
  nfirst[1] = 1;
  nlast[1] = 1;
  nat[1] = 6;
  c[1][1] = 1.0;
  nclose = 1;
  nopen = 1;
  fract = 1.0;
  p.assign(2, 0.0);
  coord[1][1] = 0.5;
  coord[2][1] = 0.5;
  coord[3][1] = 0.5;   // avoid origin / zero planes: x/r^3, y/r^3, z/r^3 columns
  tore[6] = 4.0;
  for (int ig = 1; ig <= 6; ++ig) {
    zzz[1][ig] = 0.0;
    ccc[1][ig] = 0.0;
  }
  keywrd = " CUBE ";
  jobnam = "test_c9";
  std::remove("test_c9.grd");
  new_esp();
  std::ifstream f("test_c9.grd");
  check("T5 grd file created", f.good());
  std::string line;
  std::getline(f, line);
  check("T5 header ' 4 Density'", line.find("4 Density") != std::string::npos);
  int nlines = 0;
  bool has_inf = false;
  while (std::getline(f, line)) {
    ++nlines;
    if (line.find("inf") != std::string::npos || line.find("nan") != std::string::npos)
      has_inf = true;
  }
  check("T5 grd non-empty body", nlines > 100);
  check("T5 grd values finite", !has_inf);
  check("T5 no moperr", !moperr);
}

int main() {
  setup();
  t_rys();
  setup();
  t_vint();
  setup();
  t_qrls();
  setup();
  t_density();
  setup();
  t_new_esp();
  if (failures == 0) {
    std::cout << "ALL PASS\n";
    return 0;
  }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
