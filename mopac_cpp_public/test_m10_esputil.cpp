// test_m10_esputil.cpp — batch tests for the esp_utilities + new_esp
// subsystem (M10 batch 6): rys, vint, evec, BLAS family,
// get_minus_point_five_overlap, setupg, and an end-to-end new_esp run.
#define _CRT_SECURE_NO_WARNINGS
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "chanel_C.h"
#include "esp_C.h"
#include "esp_support.h"
#include "esp_utilities.h"
#include "funcon_C.h"
#include "matout.h"
#include "molkst_C.h"
#include "new_esp.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "setupg.h"

using common_arrays_C::coord;
using common_arrays_C::c;
using common_arrays_C::p;
using common_arrays_C::nat;
using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using common_arrays_C::h;
using molkst_C::numat;
using molkst_C::norbs;
using molkst_C::nclose;
using molkst_C::nopen;
using molkst_C::fract;
using molkst_C::keywrd;
using molkst_C::jobnam;
using parameters_C::tore;
using parameters_C::zs;
using parameters_C::zp;
using parameters_C::betas;
using parameters_C::betap;
using parameters_C::betad;
using esp_C::ixn;
using esp_C::iyn;
using esp_C::izn;
using esp_C::jxn;
using esp_C::jyn;
using esp_C::jzn;
using overlaps_C::ccc;
using overlaps_C::zzz;
using overlaps_C::allz;
using overlaps_C::allc;

static int g_checks = 0;
static double g_err = 0.0;

static void chk(int cond, const char* msg, double rel = 0.0) {
  ++g_checks;
  if (!cond) {
    std::fprintf(stderr, "[FAIL] %s (rel err %.3e)\n", msg, rel);
    std::exit(1);
  }
}
static void chk_rel(double got, double want, double tol, const char* msg) {
  double rel = std::fabs(got - want) / (std::fabs(want) > 1e-30 ? std::fabs(want) : 1.0);
  g_err = rel > g_err ? rel : g_err;
  chk(rel < tol, msg, rel);
}

int main() {
  std::fprintf(stderr, "[t1] rys\n");
  {
    double u[4], w[4];
    // x=50 (>40 branch), nroots=2 — pure rational fits.
    rys(50.0, 2, u, w);
    chk_rel(u[1], 0.0055355765, 1e-5, "rys x=50 u1");
    chk_rel(u[2], 0.0576357518, 1e-5, "rys x=50 u2");
    chk_rel(w[1], 0.1138320423, 1e-5, "rys x=50 w1");
    chk_rel(w[2], 0.0114993715, 1e-5, "rys x=50 w2");
    chk_rel(w[1] + w[2], 0.1253314137, 1e-5, "rys x=50 w1+w2");
    // nroots=1.
    rys(50.0, 1, u, w);
    chk_rel(u[1], 0.5 / (50.0 - 0.5), 2e-9, "rys x=50 nroots=1 u1");
    chk_rel(w[1], 0.1253314137, 1e-5, "rys x=50 nroots=1 w1");
    // x=20 (15..33 branch) — w1 dominated by sqrt(pie4/x)=sqrt(pi/4/20).
    rys(20.0, 1, u, w);
    chk_rel(w[1], 0.198164, 2e-4, "rys x=20 w1");
    // x=0.1 (1..3e-7 branch): w1 = 2 x f1 + exp(-x), monotone positive.
    rys(0.1, 1, u, w);
    chk(w[1] > 0.0 && u[1] > 0.0, "rys x=0.1 positive");
    // x -> 0 (NOT verified branch).
    rys(1.0e-9, 1, u, w);
    chk_rel(u[1], 0.5 - 1.0e-9 / 5.0, 1e-12, "rys tiny u1");
    chk_rel(w[1], 1.0 - 1.0e-9 / 3.0, 1e-12, "rys tiny w1");
    // nroots=2 sum invariant across branches.
    rys(0.1, 2, u, w);
    chk(w[1] + w[2] > 0.0, "rys x=0.1 nroots=2 positive");
  }

  std::fprintf(stderr, "[t2] vint\n");
  {
    double xint, yint, zint;
    // s-s (ni=nj=1 -> k=l=1): px=1 -> xint = sum weights = sqrt(pi).
    vint(xint, yint, zint, 1, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
    chk_rel(xint, 1.772453850905, 1e-6, "vint s-s xint=sqrt(pi)");
    chk_rel(yint, 1.772453850905, 1e-6, "vint s-s yint=sqrt(pi)");
    chk_rel(zint, 1.772453850905, 1e-6, "vint s-s zint=sqrt(pi)");
    // p-p, k=l=2 (shell powers 1,1), all offsets 0, t=0.5:
    //   xint = int (0.5 s)^1 (0.5 s)^1 e^{-s^2} ds
    //        = 0.25 * E[s^2] = 0.25 * 0.5 * sqrt(pi) = 0.221557.
    vint(xint, yint, zint, 2, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    chk_rel(xint, 0.221557, 5e-6, "vint p-p xint");
    chk_rel(yint, 0.221557, 5e-6, "vint p-p yint");
    chk_rel(zint, 0.221557, 5e-6, "vint p-p zint");
    // Mixed: k=1, l=2 -> E[(0.5s)^1 (0.5s)^2] = 0.125 * E[s^3] = 0.
    vint(xint, yint, zint, 1, 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5);
    chk_rel(xint, 0.0, 1e-12, "vint k=1 l=2 odd moment 0");
  }

  std::fprintf(stderr, "[t3] evec\n");
  {
    std::vector<std::vector<double>> coord3(4, std::vector<double>(2, 0.0));
    coord3[1][1] = 0.0; coord3[2][1] = 0.0; coord3[3][1] = 0.0;
    std::vector<float> av(8, 0.0f);
    evec(av, 1.0, 2.0, 2.0, coord3, 1);
    // R^2 = 1+4+4 = 9; shellR2i = 1/9.0000001 ~ 0.111111.
    chk_rel(av[0], 0.3333333333, 1e-5, "evec shellRi");
    chk_rel(av[1], 1.0 * 0.037037037, 1e-4, "evec u*sR3i");
    chk_rel(av[4], 0.1111111111, 1e-5, "evec shellR2i");
    chk_rel(av[5], 0.037037037, 1e-4, "evec shellR3i");
    chk_rel(av[6], 0.012345679, 1e-4, "evec shellR2i^2");
  }

  std::fprintf(stderr, "[t4] BLAS family\n");
  {
    std::vector<float> a = {0, 1, 2, 3}, b = {0, 4, 5, 6};
    chk_rel(sdot(3, a, 1, b, 1), 32.0f, 1e-6, "sdot");
    saxpy(3, 2.0f, a, 1, b, 1);
    chk_rel(b[1], 6.0f, 1e-6, "saxpy b1"); chk_rel(b[3], 12.0f, 1e-6, "saxpy b3");
    sscal(3, 0.5f, a, 1);
    chk_rel(a[1], 0.5f, 1e-6, "sscal a1"); chk_rel(a[3], 1.5f, 1e-6, "sscal a3");
    chk_rel(snrm2(3, a, 1), std::sqrt(0.25f + 1.0f + 2.25f), 1e-5, "snrm2");
    std::vector<float> x = {0, 1, 2, 3}, y = {0, 7, 8, 9};
    sswap(3, x, 1, y, 1);
    chk_rel(x[1], 7.0f, 1e-6, "sswap x1"); chk_rel(y[3], 3.0f, 1e-6, "sswap y3");
    scopy(3, x, 1, y, 1);
    chk_rel(y[2], 8.0f, 1e-6, "scopy y2");
    // strided: dot(a[1],a[3]) with incx=2.
    std::vector<float> s = {0, 1, 2, 3, 4};
    chk_rel(sdot(2, s, 2, s, 2), 1.0f * 1.0f + 3.0f * 3.0f, 1e-6, "sdot strided");
  }

  std::fprintf(stderr, "[t5] setupg + get_minus_point_five_overlap\n");
  {
    numat = 1; norbs = 1;
    nat.assign(3, 0); nfirst.assign(3, 0); nlast.assign(3, 0);
    nat[1] = 1; nfirst[1] = 1; nlast[1] = 1;
    zs[1] = 1.2458568; zp[1] = 0.0;
    betas[1] = 1.0; betap[1] = 0.0; betad[1] = 0.0;
    setupg();
    chk_rel(ccc[1][1], 9.163596280e-3, 1e-7, "setupg ccc H 1s c1");
    chk_rel(ccc[1][6], 1.303340841e-1, 1e-7, "setupg ccc H 1s c6");
    double wantZ = 23.10303149 * zs[1] * zs[1];
    chk_rel(zzz[1][1], wantZ, 1e-6, "setupg zzz H 1s z1");
    // H: 1 orbital; h = {1} -> S^{-1/2} = 1.
    h.assign(4, 0.0); h[1] = 1.0;
    std::vector<std::vector<double>> s(2, std::vector<double>(2, 0.0));
    get_minus_point_five_overlap(s);
    chk_rel(s[1][1], 1.0, 1e-12, "S^-1/2 H s11");
    // Two-orbital check: S = [[1,0.5],[0.5,1]] -> S^{-1/2} symmetric,
    // and (S^{-1/2}) S (S^{-1/2}) = I.
    numat = 2; norbs = 2;
    nat[2] = 1; nfirst[2] = 2; nlast[2] = 2;
    h[1] = 1.0; h[2] = 0.5; h[3] = 1.0;
    std::vector<std::vector<double>> s2(3, std::vector<double>(3, 0.0));
    get_minus_point_five_overlap(s2);
    chk_rel(s2[1][2], s2[2][1], 1e-12, "S^-1/2 symmetric");
    // Rebuild S from the standard 1s-1s form used here: S = [1,0.5;0.5,1].
    double S[3][3] = {{0, 0, 0}, {0, 1.0, 0.5}, {0, 0.5, 1.0}};
    // Direct check: T = S^{-1/2} S S^{-1/2} = I
    double T[3][3] = {{0}};
    for (int a = 1; a <= 2; ++a)
      for (int k = 1; k <= 2; ++k) {
        double t1 = 0.0;
        for (int i = 1; i <= 2; ++i) t1 += s2[a][i] * S[i][k];
        for (int b = 1; b <= 2; ++b) T[a][b] += t1 * s2[k][b];
      }
    chk_rel(T[1][1], 1.0, 1e-9, "S^-1/2 S S^-1/2 = I (1,1)");
    chk_rel(T[1][2], 0.0, 1e-9, "S^-1/2 S S^-1/2 = I (1,2)");
    chk_rel(T[2][2], 1.0, 1e-9, "S^-1/2 S S^-1/2 = I (2,2)");
  }

  std::fprintf(stderr, "[t5b] matout\n");
  {
    using chanel_C::iw;
    iw = 6;
    // 3x3 column-major (0-based) eigenvectors + 1-based eigenvalues.
    double av[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double ev[4] = {0, 1.5, 2.5, 3.5};
    int nr = 3;
    std::fflush(stdout);
    FILE* cap = std::freopen("_matout_cap.txt", "w", stdout);
    if (cap) {
      matout(av, ev, 3, nr, 3);
      std::fflush(stdout);
      std::fclose(stdout);
    }
    std::string out;
    FILE* f = std::fopen("_matout_cap.txt", "r");
    if (f) {
      char buf[2048];
      while (std::fgets(buf, sizeof(buf), f)) out += buf;
      std::fclose(f);
    }
    chk(out.find("ROOT NO.") != std::string::npos, "matout ROOT NO line");
    chk(out.find("1.50000") != std::string::npos, "matout eigenvalue printed");
    chk(out.find("1.00000") != std::string::npos, "matout eigenvector printed");
    std::remove("_matout_cap.txt");
  }

  std::fprintf(stderr, "[t6] new_esp end-to-end (H atom, CUBE)\n");  {
    numat = 1; norbs = 1; nclose = 1; nopen = 1; fract = 1.0;
    nat.assign(2, 0); nfirst.assign(2, 0); nlast.assign(2, 0);
    nat[1] = 1; nfirst[1] = 1; nlast[1] = 1;
    coord.assign(4, std::vector<double>(2, 0.0));
    c.assign(2, std::vector<double>(2, 0.0));
    c[1][1] = 1.0;
    p.assign(2, 0.0);
    h.assign(2, 0.0); h[1] = 1.0;
    zs[1] = 1.2458568; zp[1] = 0.0;
    tore[1] = 1.0;
    betas[1] = 1.0; betap[1] = 0.0; betad[1] = 0.0;
    setupg();
    keywrd = " CUBE";
    jobnam = "test_m10_esputil";
    std::remove((jobnam + ".grd").c_str());
    new_esp();
    std::ifstream f(jobnam + ".grd");
    chk(f.good(), "new_esp CUBE file written");
    std::string line;
    std::getline(f, line);
    chk(line.find(" 4 Density") != std::string::npos, "CUBE header 1");
    std::getline(f, line);
    chk(line.find("Electron density") != std::string::npos, "CUBE header 2");
    f.close();
    std::remove((jobnam + ".grd").c_str());
  }

  std::fprintf(stderr, "ALL %d CHECKS PASS (max rel err %.3e)\n", g_checks, g_err);
  return 0;
}
