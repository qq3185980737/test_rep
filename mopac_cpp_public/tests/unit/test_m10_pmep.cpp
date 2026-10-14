// test_m10_pmep.cpp — batch tests for pmep.F90 translation (M10 batch 8):
// genvec (Connolly unit vectors), collis (probe collision), drepp2 (AM1
// two-electron repulsion, H-H SS/SS path), drotat (H path: enuc = tore*gam).
#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "parameters_C.h"
#include "pmep.h"
#include "rotate_C.h"

using namespace funcon_C;
using parameters_C::natorb;
using parameters_C::am;
using parameters_C::tore;
using parameters_C::dd;
using parameters_C::qq;
using parameters_C::ad;
using parameters_C::aq;
using rotate_C::ccore;

static int g_checks = 0;
static double g_err = 0.0;
static void chk(int cond, const char* msg) {
  ++g_checks;
  if (!cond) { std::fprintf(stderr, "[FAIL] %s\n", msg); std::exit(1); }
}
static void chk_rel(double got, double want, double tol, const char* msg) {
  double rel = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
  g_err = rel > g_err ? rel : g_err;
  ++g_checks;
  if (!(rel < tol)) { std::fprintf(stderr, "[FAIL] %s (rel %.3e, got %.8f want %.8f)\n", msg, rel, got, want); std::exit(1); }
}

int main() {
  std::fprintf(stderr, "[t1] genvec unit vectors over sphere\n");
  {
    std::vector<std::vector<double>> u(4, std::vector<double>(51, 0.0));
    int n = 50;
    genvec(u, n);
    chk(n > 0 && n <= 50, "genvec n reduced");
    chk_rel(u[1][1], 0.0, 1e-12, "genvec first x=0");
    chk_rel(u[2][1], 0.0, 1e-12, "genvec first y=0");
    chk_rel(u[3][1], 1.0, 1e-12, "genvec first z=1");
    for (int i = 1; i <= n; ++i) {
      double norm2 = u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i];
      chk_rel(norm2, 1.0, 1e-10, "genvec unit norm");
    }
  }

  std::fprintf(stderr, "[t2] collis probe collision\n");
  {
    std::vector<double> cw = {0.0, 1.0, 1.0, 1.0};
    std::vector<std::vector<double>> cnbr(4, std::vector<double>(2, 0.0));
    std::vector<double> rnbr(2, 0.0);
    cnbr[1][1] = 0.5; cnbr[2][1] = 1.0; cnbr[3][1] = 1.0; rnbr[1] = 1.0;
    chk(collis(cw, 1.0, cnbr, rnbr, 1, 0), "collis hit (distance < sum radii)");
    cnbr[1][1] = 5.0; cnbr[2][1] = 5.0; cnbr[3][1] = 5.0; rnbr[1] = 1.0;
    chk(!collis(cw, 1.0, cnbr, rnbr, 1, 0), "collis miss (far)");
    chk(!collis(cw, 1.0, cnbr, rnbr, 0, 0), "collis nnbr=0 -> false");
  }

  std::fprintf(stderr, "[t3] drepp2 H-H SS/SS repulsion\n");
  {
    natorb[1] = 1;
    am[1] = 1.00797;
    tore[1] = 1.0;
    std::vector<double> ri(23, 0.0);
    double core[4][3] = {};
    drepp2(1, 1.0, ri, core);
    // ri(1) = ev/sqrt((rij/a0)^2 + (0.5/am)^4); python: 13.92779342
    chk_rel(ri[1], 13.92779342, 1e-6, "drepp2 ri1");
    chk_rel(core[1][1], ri[1], 1e-12, "drepp2 core11=ri1");
    chk_rel(core[1][2], tore[1] * ri[1], 1e-12, "drepp2 core12=tore*ri1");
  }

  std::fprintf(stderr, "[t4] drotat H path enuc=tore*gam\n");
  {
    natorb[1] = 1;
    am[1] = 1.00797;
    tore[1] = 1.0;
    std::vector<double> xi = {0.0, 0.0, 0.0, 0.0};
    std::vector<double> xj = {0.0, 1.0, 0.0, 0.0};
    std::vector<double> e1b(11, 0.0);
    double enuc = 0.0;
    drotat(1, xi, xj, e1b, enuc, 1.0);
    chk_rel(enuc, tore[1] * 13.92779342, 1e-6, "drotat enuc = tore*gam");
  }

  std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
  return 0;
}
