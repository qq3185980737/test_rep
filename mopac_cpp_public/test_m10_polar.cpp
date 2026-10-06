// test_m10_polar.cpp — batch tests for polar.F90 translation (M10 batch 11):
// independent matrix/trace/parsing helpers of the polarizability driver.
#define _CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "polar.h"
#include "polar_helpers.h"

using molkst_C::norbs;
using molkst_C::numat;
using molkst_C::keywrd;
using common_arrays_C::labels;
using parameters_C::polvol;

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
  std::fprintf(stderr, "[t1] zerom\n");
  {
    std::vector<std::vector<double>> x(4, std::vector<double>(4, 7.0));
    zerom(x, 3);
    for (int i = 1; i <= 3; ++i)
      for (int j = 1; j <= 3; ++j) chk(x[i][j] == 0.0, "zerom zeroed");
    chk(x[3][3] == 0.0, "zerom 3x3");
  }

  std::fprintf(stderr, "[t2] tf (T matrix from ua/ga/ub/gb)\n");
  {
    norbs = 2;
    // logical 1-based rows live at physical index i (index 0 unused)
    std::vector<std::vector<double>> ua = {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> ga = {{0, 0, 0}, {0, 0, 1}, {0, 0, 0}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> ub = {{0, 0, 0}, {0, 0, 0}, {0, 1, 0}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> gb = {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}, {0, 0, 0, 0}};
    // python: t11=1, t12=0, t21=0, t22=-1
    std::vector<std::vector<double>> t(3, std::vector<double>(3, 0.0));
    tf(ua, ga, ub, gb, t, 2);
    chk_rel(t[1][1], 1.0, 1e-12, "tf t11=1");
    chk_rel(t[1][2], 0.0, 1e-12, "tf t12=0");
    chk_rel(t[2][1], 0.0, 1e-12, "tf t21=0");
    chk_rel(t[2][2], -1.0, 1e-12, "tf t22=-1");
  }

  std::fprintf(stderr, "[t3] transf (C^T f C)\n");
  {
    norbs = 2;
    std::vector<std::vector<double>> f = {{0, 0, 0}, {0, 1, 2}, {0, 3, 4}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> c = {{0, 0, 0}, {0, 2, 0}, {0, 0, 1}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> g(3, std::vector<double>(3, 0.0));
    transf(f, g, c, 2);
    // python: g11=4, g12=4, g21=6, g22=4
    chk_rel(g[1][1], 4.0, 1e-12, "transf g11=4");
    chk_rel(g[1][2], 4.0, 1e-12, "transf g12=4");
    chk_rel(g[2][1], 6.0, 1e-12, "transf g21=6");
    chk_rel(g[2][2], 4.0, 1e-12, "transf g22=4");
  }

  std::fprintf(stderr, "[t4] trace variants (trsub/trudgu/trugdu/trugud)\n");
  {
    std::vector<std::vector<double>> ul = {{0, 0, 0}, {0, 1, 2}, {0, 3, 4}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> x = {{0, 0, 0}, {0, 5, 6}, {0, 7, 8}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> ur = {{0, 0, 0}, {0, 2, 0}, {0, 1, 1}, {0, 0, 0, 0}};
    // python: trsub=220, trudgu=252, trugdu=220, trugud=262
    chk_rel(trsub(ul, x, ur, 2, 2, 2), 220.0, 1e-10, "trsub=220");
    chk_rel(trudgu(ul, x, ur, 2, 2, 2), 252.0, 1e-10, "trudgu=252");
    chk_rel(trugdu(ul, x, ur, 2, 2, 2), 220.0, 1e-10, "trugdu=220");
    chk_rel(trugud(ul, x, ur, 2, 2, 2), 262.0, 1e-10, "trugud=262");
  }

  std::fprintf(stderr, "[t5] wrdkey POLAR(...) keyword parsing\n");
  {
    keywrd = "POLAR(IWFLB=2,BETA=1,TOL=0.001)";
    chk_rel(wrdkey(keywrd, "POLAR(", 6, "IWFLB", 5, 0.0), 2.0, 1e-12, "wrdkey IWFLB=2");
    chk_rel(wrdkey(keywrd, "POLAR(", 6, "BETA", 4, 1.0), 1.0, 1e-12, "wrdkey BETA=1");
    chk_rel(wrdkey(keywrd, "POLAR(", 6, "TOL", 3, 0.001), 0.001, 1e-6, "wrdkey TOL");
    chk_rel(wrdkey(keywrd, "POLAR(", 6, "NOKEY", 5, 42.0), 42.0, 1e-12, "wrdkey default");
  }

  std::fprintf(stderr, "[t6] aval\n");
  {
    norbs = 2;
    std::vector<std::vector<double>> h = {{0, 0, 0}, {0, 1, 2}, {0, 3, 4}, {0, 0, 0, 0}};
    std::vector<std::vector<double>> d = {{0, 0, 0}, {0, 5, 6}, {0, 7, 8}, {0, 0, 0, 0}};
    // python: -sum(h_ij d_ji) = -69
    chk_rel(aval(h, d, 2), -69.0, 1e-12, "aval=-69");
  }

  std::fprintf(stderr, "[t7] pol_vol\n");
  {
    numat = 1;
    labels.assign(2, 0);
    labels[1] = 6;
    std::fill(polvol, polvol + 101, 0.0);
    // No polvol override -> polarizability/a0^3 = average.
    chk_rel(pol_vol(2.5), 2.5, 1e-12, "pol_vol plain");
    // With polvol[99]=2, polvol[100]=1, polvol[labels[1]]=3:
    //   polariz = avg*a0^3*2 + 1 + 3; /a0^3
    polvol[99] = 2.0; polvol[100] = 1.0; polvol[labels[1]] = 3.0;
    double want = (2.5 * std::pow(0.5291772083, 3.0) * 2.0 + 1.0 + 3.0) /
                  std::pow(0.5291772083, 3.0);
    chk_rel(pol_vol(2.5), want, 1e-9, "pol_vol with polvol table");
  }

  std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
  return 0;
}
