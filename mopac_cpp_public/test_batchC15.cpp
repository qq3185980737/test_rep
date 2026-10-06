// test_batchC15.cpp — batch tests for 批C15 (pmep).
// T1 pmep non-AM1 branch; T2 pmepco single-H numerical check;
// T3 genvec unit vectors; T4 collis; T5 drepp2 H integral;
// T6 packp; T7 grids single H; T8 surfa single H; T9 mepchg single point;
// T10 mepmap single atom completes.
#include "pmep.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "rotate_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;

static int failures = 0;

static void check(const char* name, bool ok) {
  if (!ok) {
    std::cerr << "FAIL: " << name << "\n";
    ++failures;
  } else {
    std::cout << "ok: " << name << "\n";
  }
}

static void base_setup() {
  numcal = 15;
  moperr = false;
  mozyme = false;
  keywrd = " AM1";
  natoms = 1;
  numat = 1;
  nvar = 0;
  nat.assign(2, 1);
  labels.assign(2, 1);
  nfirst.assign(2, 1);
  nlast.assign(2, 1);
  na.assign(3, 0);
  nb.assign(3, 0);
  nc.assign(3, 0);
  p.assign(12, 0.0);
  p[1] = 1.0;
  std::fill(tore, tore + 2, 0.0);
  tore[1] = 1.0;
  std::fill(am, am + 2, 0.0);
  am[1] = 0.5;
  natorb[1] = 1;
  dd[1] = 0.0;
  qq[1] = 0.0;
  ad[1] = 1.0;
  aq[1] = 1.0;
  coord.assign(4, std::vector<double>(2, 0.0));
  geo.assign(4, std::vector<double>(2, 0.0));
  itemp_1 = 50 * numat + 1000;
  mpack = 100;
  lm61 = 50 * numat + 1000;
  ch.assign(2, 0.0);
}

int main() {
  std::cout << "T1 start" << std::endl;
  // T1: pmep without AM1 -> "PMEP WAS NOT EXECUTED", returns quickly
  base_setup();
  keywrd = " PM3";
  pmep();
  check("T1 pmep non-AM1 returns", !moperr);

  std::cout << "T2 start" << std::endl;
  // T2: pmepco single H at r=1.0 -> ui = ri(1)*fnear*ev
  base_setup();
  {
    std::vector<double> pp(2, 0.0), ria(2, 0.0), w(4, 0.0);
    pp[1] = 1.0;
    ria[1] = 1.0;
    w[1] = 1.0;
    double ui = 0.0;
    pmepco(pp, ria, w, ui, coord, 1, 1);
    double fnear = 1.0 + std::exp(-(2.92 * (1.0 - 0.05)));
    double r = 1.0 / a0;
    double ri1 = ev / std::sqrt(r * r + (0.5 / 0.5) * (0.5 / 0.5));
    double expect = ri1 * fnear * ev;
    check("T2 pmepco H potential", std::fabs(ui - expect) < 1e-6 * std::fabs(expect));
  }

  std::cout << "T3 start" << std::endl;
  // T3: genvec produces unit vectors
  base_setup();
  {
    std::vector<std::vector<double>> u(4, std::vector<double>(200, 0.0));
    int n = 100;
    genvec(u, n);
    check("T3 genvec count>0", n > 0 && n <= 100);
    bool unit = true;
    for (int i = 1; i <= n; ++i) {
      double norm = std::sqrt(u[1][i] * u[1][i] + u[2][i] * u[2][i] + u[3][i] * u[3][i]);
      if (std::fabs(norm - 1.0) > 1e-6) { unit = false; break; }
    }
    check("T3 genvec unit norm", unit);
  }

  std::cout << "T4 start" << std::endl;
  // T4: collis detects overlap
  base_setup();
  {
    std::vector<double> cw(4, 0.0);
    std::vector<std::vector<double>> cnbr(4, std::vector<double>(201, 0.0));
    std::vector<double> rnbr(201, 0.0);
    cw[1] = 0.0; cw[2] = 0.0; cw[3] = 0.0;
    cnbr[1][1] = 0.5; cnbr[2][1] = 0.0; cnbr[3][1] = 0.0;
    rnbr[1] = 1.0;
    check("T4 collis true", collis(cw, 0.2, cnbr, rnbr, 1, 1));
    cnbr[1][1] = 5.0;
    check("T4 collis false", !collis(cw, 0.2, cnbr, rnbr, 1, 1));
  }

  std::cout << "T5 start" << std::endl;
  // T5: drepp2 hydrogen SS integral ri(1) = ev/sqrt(r^2 + (pp/am)^4)
  base_setup();
  {
    std::vector<double> ri(23, 0.0);
    double core[4][3] = {};
    drepp2(1, 1.0, ri, core);
    double r = 1.0 / a0;
    double expect = ev / std::sqrt(r * r + (0.5 / 0.5) * (0.5 / 0.5));
    check("T5 drepp2 H ri(1)", std::fabs(ri[1] - expect) < 1e-9 * std::fabs(expect));
    check("T5 drepp2 core(1,1)", std::fabs(core[1][1] - ri[1]) < 1e-12);
    check("T5 drepp2 core(1,2)=tore*ri", std::fabs(core[1][2] - tore[1] * ri[1]) < 1e-12);
  }

  std::cout << "T6 start" << std::endl;
  // T6: packp rewrites density matrix
  base_setup();
  {
    nfirst[1] = 1;
    nlast[1] = 2;
    p.assign(8, 0.0);
    p[1] = 0.1; p[2] = 0.2; p[3] = 0.3;
    std::vector<double> pp(12, 0.0);
    int mn = 0;
    packp(p, pp, mn);
    check("T6 packp mn=3", mn == 3);
    check("T6 packp pp[1]=p[1]", std::fabs(pp[1] - 0.1) < 1e-12);
    check("T6 packp pp[2..3]=p[2..3]", std::fabs(pp[2] - 0.2) < 1e-12 && std::fabs(pp[3] - 0.3) < 1e-12);
  }

  std::cout << "T7 start" << std::endl;
  // T7: grids single H at origin generates Williams points
  base_setup();
  {
    std::vector<std::vector<double>> potpt(4, std::vector<double>(itemp_1 + 1, 0.0));
    int nmep = 0;
    grids(coord, potpt, nmep);
    check("T7 grids nmep>0", nmep > 0);
  }

  std::cout << "T8 start" << std::endl;
  // T8: surfa single H generates Connolly points
  base_setup();
  {
    std::vector<std::vector<double>> potpt(4, std::vector<double>(itemp_1 + 1, 0.0));
    int nmep = 0;
    surfa(coord, potpt, nmep);
    check("T8 surfa nmep>0", nmep > 0);
  }

  std::cout << "T9 start" << std::endl;
  // T9: mepchg with a single sample point completes and prints
  base_setup();
  {
    std::vector<std::vector<double>> potpt(4, std::vector<double>(2, 0.0));
    potpt[1][1] = 2.0; potpt[2][1] = 0.0; potpt[3][1] = 0.0;
    mepchg(coord, potpt, 1);
    check("T9 mepchg completes", !moperr);
  }

  std::cout << "T10 start" << std::endl;
  // T10: mepmap single atom completes (no minimum expected for repulsive H)
  base_setup();
  {
    std::vector<std::vector<double>> c1(4, std::vector<double>(2, 0.0));
    mepmap(c1);
    check("T10 mepmap completes", !moperr);
  }

  if (failures == 0) {
    std::cout << "ALL PASS\n";
    return 0;
  }
  std::cerr << failures << " FAILURES\n";
  return 1;
}
