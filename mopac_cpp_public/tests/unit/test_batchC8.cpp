// test_batchC8.cpp — regression tests for react1 + dock (batch C8).
#include "react1.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;

static int cal_ctr = 0;

static void setup_geom() {
  numcal = ++cal_ctr;
  keywrd = " ";
  moperr = false;
  mozyme = false;
  numat = 3;
  natoms = 3;
  id = 0;
  nvar = 9;
  geo.assign(4, std::vector<double>(4, 0.0));
  geoa.assign(4, std::vector<double>(4, 0.0));
  coord.assign(4, std::vector<double>(4, 0.0));
  labels.assign(4, 6);
  nat.assign(4, 6);
  na.assign(4, 0);
  nb.assign(4, 0);
  nc.assign(4, 0);
  xparam.assign(20, 0.0);
  grad.assign(20, 0.0);
  p.assign(20, 0.0);
  pa.assign(20, 0.0);
  pb.assign(20, 0.0);
  step = 0.0;
  tleft = 0.0;
  escf = 0.0;
  gnorm = 1.0;
  cosine = 1.0;
  iflepo = 0;
}

static std::string read_file(const char* path) {
  FILE* f = std::fopen(path, "rb");
  if (!f) return "";
  std::string s;
  char buf[4096];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) s.append(buf, n);
  std::fclose(f);
  return s;
}

static void T1_dock_minimizes() {
  // geo: triangle with centroid at origin (so dock's saved sum_a == 0).
  setup_geom();
  geo[1][1] = 1.0;  geo[2][1] = 1.0;  geo[3][1] = 0.0;
  geo[1][2] = -1.0; geo[2][2] = 0.0;  geo[3][2] = 0.0;
  geo[1][3] = 0.0;  geo[2][3] = -1.0; geo[3][3] = 0.0;
  // geoa: same shape, rotated 30 deg about z, then translated.
  const double ca = std::cos(30.0 * 3.141592653589793 / 180.0);
  const double sa = std::sin(30.0 * 3.141592653589793 / 180.0);
  for (int i = 1; i <= 3; ++i) {
    double x = geo[1][i], y = geo[2][i];
    geoa[1][i] = ca * x - sa * y + 1.0;
    geoa[2][i] = sa * x + ca * y + 2.0;
    geoa[3][i] = geo[3][i] + 3.0;
  }
  double dist = -1.0;
  dock(geoa, geo, dist);
  if (dist < 0.0 || dist > 0.05) {
    std::fprintf(stderr, "FAIL T1 dock dist=%f expect small\n", dist);
    std::exit(1);
  }
}

static void T2_react1_identical_geometries() {
  // GEO_REF keyword + identical geometries -> early mopend.
  setup_geom();
  keywrd = " GEO_REF ";
  geo[1][1] = 0.0; geo[2][1] = 0.0; geo[3][1] = 0.0;
  geo[1][2] = 1.0; geo[2][2] = 0.0; geo[3][2] = 0.0;
  geo[1][3] = 0.0; geo[2][3] = 1.0; geo[3][3] = 0.0;
  geoa = geo;
  coord = geo;

  if (!std::freopen("test_c8_t2.txt", "w", stdout)) {
    std::fprintf(stderr, "FAIL cannot redirect stdout\n"); std::exit(1);
  }
  react1();
  std::fflush(stdout);
  std::string out = read_file("test_c8_t2.txt");
  if (out.find("CARTESIAN GEOMETRY OF FIRST SYSTEM") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing first-system header\n"); std::exit(1);
  }
  if (out.find("THE TWO GEOMETRIES ARE IDENTICAL OR ALMOST IDENTICAL") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing identical-geometry message\n"); std::exit(1);
  }
  if (out.find("THE TWO GEOMETRIES IN SADDLE ARE IDENTICAL.") == std::string::npos) {
    std::fprintf(stderr, "FAIL T2 missing mopend message\n"); std::exit(1);
  }
  if (!moperr) {
    std::fprintf(stderr, "FAIL T2 moperr not set\n"); std::exit(1);
  }
}

int main() {
  T1_dock_minimizes();
  T2_react1_identical_geometries();
  std::fprintf(stderr, "ALL PASS\n");
  return 0;
}
