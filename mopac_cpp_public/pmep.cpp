// pmep.cpp — C++ translation of "pmep.F90" (12 program units, 1600 lines).
#define _CRT_SECURE_NO_WARNINGS
//
// PMEP — Parametric Molecular Electrostatic Potential and MEP charges (AM1).
// pmepco computes the potential at arbitrary points; surfa/grids generate
// Connolly/Williams sample surfaces; mepchg fits MEP charges; mepmap maps
// minima and contour lines; meprot selects the mapping plane; packp rewrites
// the density matrix.

#include "pmep.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "blas1.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "gmetry.h"
#include "molkst_C.h"
#include "mopend.h"
#include "osinv.h"
#include "parameters_C.h"
#include "reada.h"
#include "rotate_C.h"
#include "second.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;
using namespace rotate_C;

namespace {

// AM1 screening exponents alpw(1..107) and aw(1..107) from the F90 DATA stmts.
const double alpw_table[108] = {
    0.0, 2.92, 0.0, 0.0, 0.0, 0.0, 4.36, 6.95, 9.22, 5.53,
    0.0, 0.0,  0.0, 0.0, 0.0, 2.3, 2.6,  3.32, 0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 3.32, 0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0,
};
const double aw_table[108] = {
    0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.63, 0.64, 0.67, 0.29,
    0.0, 0.0,  0.0, 0.0, 0.45, 0.37, 0.31, 0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.31, 0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, 0.0,  0.0, 0.0, 0.0,  0.0,  0.0,
};

std::string upcase(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

}  // namespace

void pmep() {
  static int nmep = 0;
  // CHECK FOR SAVING TIME
  if (keywrd.find(" AM1") == std::string::npos) {
    std::printf(" AM1 PARAMETERS FOR H,C,N,F AND CL ARE AVAILABLE.  PMEP WAS NOT EXECUTED\n");
    return;
  }
  for (int i = 1; i <= numat; ++i)
    if (alpw_table[nat[i]] < 0.1) {
      std::printf(" AM1 PARAMETERS FOR H,C,N,F AND CL ARE AVAILABLE.  PMEP WAS NOT EXECUTED\n");
      return;
    }
  double t1 = second(1);
  gmetry(geo, coord);
  if (moperr) return;
  // CALCULATE MEP ONLY
  if (keywrd.find(" PMEP") != std::string::npos) {
    FILE* fp = fopen(mep_fn.empty() ? "mep.dat" : mep_fn.c_str(), "a");
    if (fp) fclose(fp);  // open/append semantics preserved; mepmap writes below
    mepmap(coord);
    if (moperr) return;
    t1 = second(1) - t1;
    std::printf("\n TIME TO CALCULATE MEP%9.2f SECONDS\n", t1);
    return;
  }
  itemp_1 = 50 * numat + 1000;
  // CHOOSE SURFACE WHERE POINTS ARE SAMPLED
  std::vector<std::vector<double>> potpt(4, std::vector<double>(itemp_1 + 1, 0.0));
  if (keywrd.find("WILLIAMS") != std::string::npos) {
    // THE WILLIAMS SURFACE
    grids(coord, potpt, nmep);
    if (moperr) return;
  } else {
    // THE CONNOLLY SURFACE
    surfa(coord, potpt, nmep);
    if (moperr) return;
  }
  // CALCULATE MEP CHARGES
  mepchg(coord, potpt, nmep);
  if (moperr) return;
  t1 = second(1) - t1;
  std::printf("\n TIME TO CALCULATE MEP%9.2f SECONDS\n", t1);
}

void pmepco(const std::vector<double>& pp, const std::vector<double>& ria,
            const std::vector<double>& w, double& ui,
            const std::vector<std::vector<double>>& co, int nonzo, int id) {
  std::vector<double> e1b(11, 0.0);
  std::vector<double> vp(static_cast<std::size_t>(nonzo) + 1, 0.0);
  double enuclr = 0.0;
  int mn = 0;
  for (int i = 1; i <= numat; ++i) {
    int ia = nfirst[i];
    int ib = nlast[i];
    int ni = nat[i];
    double rr;
    if (id == 1) {
      rr = ria[i];
    } else {
      double xx = w[1] - co[1][i];
      double yy = w[2] - co[2][i];
      double zz = w[3] - co[3][i];
      rr = std::sqrt(xx * xx + yy * yy + zz * zz);
    }
    double fnear = 1.0 + std::exp(-(alpw_table[ni] * (rr - aw_table[ni])));
    double enuc = 0.0;
    // Fortran passes co(1,i): the three coordinates of atom i (column start).
    std::vector<double> xi4{0.0, co[1][i], co[2][i], co[3][i]};
    drotat(ni, xi4, w, e1b, enuc, rr);
    enuclr = enuclr + enuc * fnear;
    int i2 = 0;
    for (int i1 = ia; i1 <= ib; ++i1) {
      if (i1 - ia + 1 > 0) {
        for (int t = 0; t < i1 - ia + 1; ++t) vp[mn + 1 + t] = e1b[i2 + 1 + t];
        i2 = i1 - ia + 1 + i2;
        mn = i1 - ia + 1 + mn;
      }
    }
  }
  mn = 0;
  ui = 0.0;
  for (int i = 1; i <= numat; ++i) {
    int ia = nfirst[i];
    int ib = nlast[i];
    for (int j = ia; j <= ib; ++j) {
      for (int k = ia; k <= j; ++k) {
        mn = mn + 1;
        if (k == j) {
          ui = ui + pp[mn] * vp[mn];
        } else {
          ui = ui + 2.0 * pp[mn] * vp[mn];
        }
      }
    }
  }
  ui = (enuclr + ui) * ev;
  // UI IN KCAL/MOL
}

void surfa(const std::vector<std::vector<double>>& co,
           std::vector<std::vector<double>>& potpt, int& nmep) {
  double vander[101] = {};
  vander[1] = 1.20; vander[2] = 1.20; vander[3] = 1.37; vander[4] = 1.45;
  vander[5] = 1.45; vander[6] = 1.50; vander[7] = 1.50; vander[8] = 1.40;
  vander[9] = 1.35; vander[10] = 1.30; vander[11] = 1.57; vander[12] = 1.36;
  vander[13] = 1.24; vander[14] = 1.17; vander[15] = 1.80; vander[16] = 1.75;
  vander[17] = 1.70; vander[35] = 2.3;
  double scale = 1.4, dens = 1.0, scincr = 0.20;
  int layer = 4;
  itemp_1 = 50 * numat + 1000;
  int mxmep = itemp_1;
  std::vector<int> ias(mxmep + 1, 0);
  std::vector<double> rad(mxmep + 1, 0.0);
  double pi = 4.0 * std::atan(1.0);
  vander[30] = 1.00;  // INSERT VAN DER WAAL RADII FOR ZINC

  if (keywrd.find("SCALE=") != std::string::npos)
    scale = reada(keywrd, static_cast<int>(keywrd.find("SCALE=")));
  if (keywrd.find("DEN=") != std::string::npos)
    dens = reada(keywrd, static_cast<int>(keywrd.find("DEN=")));
  if (keywrd.find("SCINCR=") != std::string::npos)
    scincr = reada(keywrd, static_cast<int>(keywrd.find("SCINCR=")));
  if (keywrd.find("NSURF=") != std::string::npos)
    layer = static_cast<int>(std::lround(reada(keywrd, static_cast<int>(keywrd.find("NSURF=")))));

  std::vector<double> ci(4, 0.0), temp0(4, 0.0);
  std::vector<std::vector<double>> cw(4, std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> cnbr(4, std::vector<double>(201, 0.0));
  std::vector<double> rnbr(201, 0.0);
  std::vector<std::vector<double>> con(4, std::vector<double>(1001, 0.0));

  for (int ishel = 1; ishel <= layer; ++ishel) {
    double rw = 0.0;
    double den = dens;
    for (int i = 1; i <= numat; ++i) {
      int ipoint = nat[i];
      rad[i] = vander[ipoint] * scale;
      if (rad[i] < 0.01)
        std::printf(" VAN DER WAALS' RADIUS FOR ATOM %3d IS ZERO, SUPPLY A VALUE IN SUBROUTINE SURFAC)\n", i);
      ias[i] = 2;
    }
    // BIG LOOP FOR EACH ATOM
    for (int iatom = 1; iatom <= numat; ++iatom) {
      if (ias[iatom] == 0) continue;
      double ri = rad[iatom];
      bool si = ias[iatom] == 2;
      for (int d = 1; d <= 3; ++d) ci[d] = co[d][iatom];
      // GATHER THE NEIGHBORING ATOMS OF IATOM
      int nnbr = 0;
      for (int jatom = 1; jatom <= numat; ++jatom) {
        if (iatom == jatom || ias[jatom] == 0) continue;
        double d2 = (ci[1] - co[1][jatom]) * (ci[1] - co[1][jatom]) +
                    (ci[2] - co[2][jatom]) * (ci[2] - co[2][jatom]) +
                    (ci[3] - co[3][jatom]) * (ci[3] - co[3][jatom]);
        if (d2 >= (2 * rw + ri + rad[jatom]) * (2 * rw + ri + rad[jatom])) continue;
        nnbr = nnbr + 1;
        if (nnbr > 200) {
          std::printf("ERROR%2sTOO MANY NEIGHBORS:%5d\n", "", nnbr);
          mopend("TOO MANY NEIGHBORS");
          return;
        }
        for (int d = 1; d <= 3; ++d) cnbr[d][nnbr] = co[d][jatom];
        rnbr[nnbr] = rad[jatom];
      }
      // CONTACT SURFACE
      if (!si) continue;
      int ncon = static_cast<int>((4 * pi * ri * ri) * den);
      ncon = std::min(1000, ncon);
      if (ncon == 0) {
        std::printf(" VECTOR LENGTH OF ZERO IN SURFAC\n");
        mopend("VECTOR LENGTH OF ZERO IN SURFAC");
        return;
      }
      genvec(con, ncon);
      // CONTACT PROBE PLACEMENT LOOP
      for (int i = 1; i <= ncon; ++i) {
        for (int d = 1; d <= 3; ++d) cw[d][1] = ci[d] + (ri + rw) * con[d][i];
        // CHECK FOR COLLISION WITH NEIGHBORING ATOMS
        std::vector<double> cw4{0.0, cw[1][1], cw[2][1], cw[3][1]};
        if (collis(cw4, rw, cnbr, rnbr, nnbr, 1)) continue;
        for (int d = 1; d <= 3; ++d) temp0[d] = ci[d] + ri * con[d][i];
        // STORE POINT IN POTPT AND INCREMENT NMEP
        nmep = nmep + 1;
        if (nmep > mxmep) goto l100;
        potpt[1][nmep] = temp0[1];
        potpt[2][nmep] = temp0[2];
        potpt[3][nmep] = temp0[3];
      }
    }
    scale = scale + scincr;
  }
  return;
l100:
  std::printf("\nERROR - TO MANY POINTS GENERATED IN SURFAC\n"
              "    REDUCE NSURF, SCALE, DEN, OR SCINCR\n");
  mopend("ERROR - TOO MANY POINTS GENERATED IN SURFAC.");
}

bool collis(const std::vector<double>& cw, double rw,
            const std::vector<std::vector<double>>& cnbr,
            const std::vector<double>& rnbr, int nnbr, int ishape) {
  if (nnbr <= 0) return false;
  if (ishape == 3) {
  } else {
    for (int i = 1; i <= nnbr; ++i) {
      double sumrad = rw + rnbr[i];
      double vect1 = std::fabs(cw[1] - cnbr[1][i]);
      if (vect1 >= sumrad) continue;
      double vect2 = std::fabs(cw[2] - cnbr[2][i]);
      if (vect2 >= sumrad) continue;
      double vect3 = std::fabs(cw[3] - cnbr[3][i]);
      if (vect3 >= sumrad) continue;
      double sr2 = sumrad * sumrad;
      double dd2 = vect1 * vect1 + vect2 * vect2 + vect3 * vect3;
      if (dd2 < sr2) return true;
    }
  }
  return false;
}

void drepp2(int ni, double rij, std::vector<double>& ri, double core[4][3]) {
  double ev1 = ev * 0.5;
  double ev2 = ev1 * 0.5;
  const double pp = 0.5;
  double r = rij / a0;
  if (natorb[ni] < 3) {
    // HYDROGEN - HYDROGEN (SS/SS)
    double aee = pp / am[ni];
    aee = aee * aee;
    ri[1] = ev / std::sqrt(r * r + aee);
    core[1][1] = ri[1];
    core[1][2] = tore[ni] * ri[1];
  } else {
    // HEAVY ATOM - HYDROGEN
    double aee = pp / am[ni];
    aee = aee * aee;
    double da = dd[ni];
    double qa = qq[ni] * 2.0;
    double ade = pp / ad[ni];
    ade = ade * ade;
    double aqe = pp / aq[ni];
    aqe = aqe * aqe;
    double rsq = r * r;
    double arg[8], sqr[8];
    arg[1] = rsq + aee;
    double xxx = r + da;
    arg[2] = xxx * xxx + ade;
    xxx = r - da;
    arg[3] = xxx * xxx + ade;
    xxx = r + qa;
    arg[4] = xxx * xxx + aqe;
    xxx = r - qa;
    arg[5] = xxx * xxx + aqe;
    arg[6] = rsq + aqe;
    arg[7] = arg[6] + qa * qa;
    for (int i = 1; i <= 7; ++i) sqr[i] = std::sqrt(arg[i]);
    double ee = ev / sqr[1];
    ri[1] = ee;
    ri[2] = ev1 / sqr[2] - ev1 / sqr[3];
    ri[3] = ee + ev2 / sqr[4] + ev2 / sqr[5] - ev1 / sqr[6];
    ri[4] = ee + ev1 / sqr[7] - ev1 / sqr[6];
    core[1][1] = ri[1];
    core[2][1] = ri[2];
    core[3][1] = ri[3];
    core[4][1] = ri[4];
    core[1][2] = tore[ni] * ri[1];
  }
}

void drotat(int ni, const std::vector<double>& xi, const std::vector<double>& xj,
            std::vector<double>& e1b, double& enuc, double rij) {
  std::vector<double> x(4, 0.0), y(4, 0.0), z(4, 0.0);
  std::vector<double> ri(23, 0.0);
  x[1] = xi[1] - xj[1];
  x[2] = xi[2] - xj[2];
  x[3] = xi[3] - xj[3];
  double rijx = rij;
  // COMPUTE INTEGRALS IN DIATOMIC FRAME
  drepp2(ni, rij, ri, ccore);
  double gam = ri[1];
  double a = 1.0 / rijx;
  x[1] = x[1] * a;
  x[2] = x[2] * a;
  x[3] = x[3] * a;
  if (std::fabs(x[3]) > 0.99999) {
    x[3] = x[3] > 0 ? 1.0 : -1.0;
    y[1] = 0.0; y[2] = 1.0; y[3] = 0.0;
    z[1] = 1.0; z[2] = 0.0; z[3] = 0.0;
  } else {
    z[3] = std::sqrt(1.0 - x[3] * x[3]);
    a = 1.0 / z[3];
    y[1] = -a * x[2] * (x[1] > 0 ? 1.0 : -1.0);
    y[2] = std::fabs(a * x[1]);
    y[3] = 0.0;
    z[1] = -a * x[1] * x[3];
    z[2] = -a * x[2] * x[3];
  }
  if (natorb[ni] > 1) {
    double xx11 = x[1] * x[1], xx21 = x[2] * x[1], xx22 = x[2] * x[2];
    double xx31 = x[3] * x[1], xx32 = x[3] * x[2], xx33 = x[3] * x[3];
    double yy11 = y[1] * y[1], yy21 = y[2] * y[1], yy22 = y[2] * y[2];
    double zz11 = z[1] * z[1], zz21 = z[2] * z[1], zz22 = z[2] * z[2];
    double zz31 = z[3] * z[1], zz32 = z[3] * z[2], zz33 = z[3] * z[3];
    double yyzz11 = yy11 + zz11;
    double yyzz21 = yy21 + zz21;
    double yyzz22 = yy22 + zz22;
    // ROTATE THE NUCLEAR ATTRACTION INTEGRALS
    e1b[1] = -css1;
    if (natorb[ni] == 4) {
      e1b[2] = -csp1 * x[1];
      e1b[3] = (-cpps1 * xx11) - cppp1 * yyzz11;
      e1b[4] = -csp1 * x[2];
      e1b[5] = (-cpps1 * xx21) - cppp1 * yyzz21;
      e1b[6] = (-cpps1 * xx22) - cppp1 * yyzz22;
      e1b[7] = -csp1 * x[3];
      e1b[8] = (-cpps1 * xx31) - cppp1 * zz31;
      e1b[9] = (-cpps1 * xx32) - cppp1 * zz32;
      e1b[10] = (-cpps1 * xx33) - cppp1 * zz33;
    }
  }
  enuc = tore[ni] * gam;
}

void genvec(std::vector<std::vector<double>>& u, int& n) {
  double pi = 4.0 * std::atan(1.0);
  int nequat = static_cast<int>(std::sqrt(n * pi));
  int nvert = nequat / 2;
  int nu = 0;
  for (int i = 1; i <= nvert + 1; ++i) {
    double fi = (pi * (i - 1)) / nvert;
    double z = std::cos(fi);
    double xy = std::sin(fi);
    int nhor = static_cast<int>(nequat * xy);
    nhor = std::max(1, nhor);
    for (int j = 1; j <= nhor; ++j) {
      double fj = (2.0 * pi * (j - 1)) / nhor;
      double x = std::cos(fj) * xy;
      double y = std::sin(fj) * xy;
      if (nu >= n) return;
      nu = nu + 1;
      u[1][nu] = x;
      u[2][nu] = y;
      u[3][nu] = z;
    }
  }
  n = nu;
}

void grids(const std::vector<std::vector<double>>& co,
           std::vector<std::vector<double>>& potpt, int& nmep) {
  double vderw[54] = {};
  vderw[1] = 2.4; vderw[5] = 3.0; vderw[6] = 2.9; vderw[7] = 2.7;
  vderw[8] = 2.6; vderw[9] = 2.55; vderw[15] = 3.1; vderw[16] = 3.05;
  vderw[17] = 3.0; vderw[35] = 3.15; vderw[53] = 3.35;
  double shell = 1.2, grid = 0.8, closer = 0.0;
  // CHECK IF VDERW IS DEFINED FOR ALL ATOMS
  int i = 0;
  for (i = 1; i <= numat; ++i)
    if (vderw[nat[i]] == 0.0) goto l20;
  goto l30;
l20:
  std::printf(" VAN DER WAALS' RADIUS NOT DEFINED FOR ATOM%d\n", i);
  std::printf(" IN WILLIAMS SURFACE ROUTINE PDGRID!\n");
  mopend("VAN DER WAALS' RADIUS NOT DEFINED IN WILLIAMS SURFACE ROUTINE PDGRID!.");
  return;
  // NOW CREATE LIMITS FOR A BOX
l30: {
  double xmin[4] = {0, 1e5, 1e5, 1e5}, xmax[4] = {0, -1e5, -1e5, -1e5};
  for (int ia = 1; ia <= numat; ++ia) {
    for (int d = 1; d <= 3; ++d) {
      if (co[d][ia] - xmin[d] < 0.0) xmin[d] = co[d][ia];
      if (co[d][ia] - xmax[d] > 0.0) xmax[d] = co[d][ia];
    }
  }
  // ADD (OR SUBTRACT) THE MAXIMUM VDERW PLUS SHELL
  double vdmax = 0.0;
  for (int d = 1; d <= 53; ++d) vdmax = std::max(vderw[d], vdmax);
  for (int d = 1; d <= 3; ++d) {
    xmin[d] = xmin[d] - vdmax - shell;
    xmax[d] = xmax[d] + vdmax + shell;
  }
  // STEP GRID BACK FROM ZERO TO FIND STARTING POINTS
  double xstart = -grid;
  while (xstart > xmin[1]) xstart = xstart - grid;
  double ystart = -grid;
  while (ystart > xmin[2]) ystart = ystart - grid;
  double zstart = -grid;
  while (zstart > xmin[3]) zstart = zstart - grid;
  std::vector<double> dist(101, 0.0);
  double zgrid = zstart;
  while (true) {
    double ygrid = ystart;
    while (true) {
      double xgrid = xstart;
      while (true) {
        bool rejected = false, accepted = false;
        for (int l = 1; l <= numat; ++l) {
          int jz = nat[l];
          dist[l] = std::sqrt((co[1][l] - xgrid) * (co[1][l] - xgrid) +
                              (co[2][l] - ygrid) * (co[2][l] - ygrid) +
                              (co[3][l] - zgrid) * (co[3][l] - zgrid));
          // REJECT GRID POINT IF ANY ATOM IS TOO CLOSE
          if (dist[l] < vderw[jz] - closer) { rejected = true; break; }
        }
        if (!rejected) {
          // BUT AT LEAST ONE ATOM MUST BE CLOSE ENOUGH
          for (int l = 1; l <= numat; ++l) {
            int jz = nat[l];
            if (dist[l] > vderw[jz] + shell) continue;
            accepted = true;
            break;
          }
        }
        if (accepted) {
          nmep = nmep + 1;
          potpt[1][nmep] = xgrid;
          potpt[2][nmep] = ygrid;
          potpt[3][nmep] = zgrid;
        }
        xgrid = xgrid + grid;
        if (xgrid <= xmax[1]) continue;
        break;
      }
      ygrid = ygrid + grid;
      if (ygrid <= xmax[2]) continue;
      break;
    }
    zgrid = zgrid + grid;
    if (zgrid <= xmax[3]) continue;
    break;
  }
}
}

void mepchg(const std::vector<std::vector<double>>& co,
            const std::vector<std::vector<double>>& potpt, int nmep) {
  static int iz = 0;
  static double rms = 0.0, rrms = 0.0, dipx = 0.0, dipy = 0.0, dipz = 0.0;
  static bool idip = false;
  double au = ev * fpc_9;
  double cf = fpc_1 * a0 * fpc_8 * 1.0e-10;
  double bohr1 = 1.0 / a0;
  double au1 = 1.0 / au;
  int numat0 = numat + 1;

  std::printf("\n\n ATOMIC CHARGES DERIVED FROM MOLECULAR ELECTROSTATIC POTENTIAL\n");
  if (keywrd.find(" WILLIAMS") != std::string::npos)
    std::printf("        (WILLIAMS SURFACES)\n\n");
  else
    std::printf("        (CONNOLLY SURFACES)\n\n");

  if (keywrd.find(" CHARGE=") != std::string::npos)
    iz = static_cast<int>(std::lround(reada(keywrd, static_cast<int>(keywrd.find(" CHARGE=")))));
  // DIPOLAR CONSTRAINTS
  if (keywrd.find("DIPOLE") != std::string::npos) {
    if (iz == 0) {
      idip = true;
      numat0 = numat + 4;
      std::printf("\n          DIPOLE CONSTRAINTS WILL BE USED\n\n");
      double dx = ux, dy = uy, dz = uz;
      if (keywrd.find("DIPX=") != std::string::npos) dx = reada(keywrd, static_cast<int>(keywrd.find("DIPX=")));
      if (keywrd.find("DIPY=") != std::string::npos) dy = reada(keywrd, static_cast<int>(keywrd.find("DIPY=")));
      if (keywrd.find("DIPZ=") != std::string::npos) dz = reada(keywrd, static_cast<int>(keywrd.find("DIPZ=")));
      dipx = dx; dipy = dy; dipz = dz;  // carry dx/dy/dz via statics
    } else {
      std::printf(" DIPOLE CONSTRAINTS NOT USED FOR CHARGED MOLECULE\n");
    }
  }
  std::vector<std::vector<double>> d(numat + 1, std::vector<double>(nmep + 1, 0.0));
  std::vector<std::vector<double>> d1(numat + 1, std::vector<double>(nmep + 1, 0.0));
  for (int k = 1; k <= numat; ++k) {
    double xk = co[1][k], yk = co[2][k], zk = co[3][k];
    for (int i = 1; i <= nmep; ++i) {
      double xp = potpt[1][i] - xk;
      double yp = potpt[2][i] - yk;
      double zp = potpt[3][i] - zk;
      double dki = std::sqrt(xp * xp + yp * yp + zp * zp);
      d[k][i] = dki;
      d1[k][i] = a0 / dki;
    }
  }

  // SET UP THE LINEAR EQUATION A*Q=B
  // MEP AT SAMPLE POINTS & B MATRIX
  int nonzo = 0;
  std::vector<double> pp(mpack + 1, 0.0);
  packp(p, pp, nonzo);
  std::vector<double> b(numat + 5, 0.0);
  std::vector<double> ep(static_cast<std::size_t>(itemp_1) + 1, 0.0);
  for (int ip = 1; ip <= nmep; ++ip) {
    double ui = 0.0;
    // Fortran passes POTPT(1,IP): the three coordinates of sample point IP.
    std::vector<double> w4{0.0, potpt[1][ip], potpt[2][ip], potpt[3][ip]};
    pmepco(pp, d[ip], w4, ui, co, nonzo, 1);
    ep[ip] = ui * au1;
    for (int i = 1; i <= numat; ++i) b[i] = b[i] + ep[ip] * d1[i][ip];
  }
  b[numat + 1] = static_cast<double>(iz);
  if (idip) {
    b[numat + 2] = dipx / cf;
    b[numat + 3] = dipy / cf;
    b[numat + 4] = dipz / cf;
  }
  // THE A(J,K) ARRAY
  std::vector<std::vector<double>> a(numat + 5, std::vector<double>(numat + 5, 0.0));
  for (int k = 1; k <= numat; ++k) {
    for (int j = 1; j <= k; ++j) {
      double ajk = 0.0;
      for (int ip = 1; ip <= nmep; ++ip) ajk += d1[k][ip] * d1[j][ip];
      a[k][j] = ajk;
      a[j][k] = ajk;
    }
    a[numat + 1][k] = 1.0;
    a[k][numat + 1] = 1.0;
    if (!idip) continue;
    a[numat + 2][k] = co[1][k] * bohr1;
    a[numat + 3][k] = co[2][k] * bohr1;
    a[numat + 4][k] = co[3][k] * bohr1;
    a[k][numat + 2] = a[numat + 2][k];
    a[k][numat + 3] = a[numat + 3][k];
    a[k][numat + 4] = a[numat + 4][k];
  }
  a[numat + 1][numat + 1] = 0.0;
  if (idip) {
    a[numat + 2][numat + 2] = 0.0;
    a[numat + 3][numat + 3] = 0.0;
    a[numat + 4][numat + 4] = 0.0;
  }
  // SOLVE AQ=B
  std::vector<double> al(numat0 * numat0, 0.0);
  int l = 0;
  for (int i = 1; i <= numat0; ++i)
    for (int j = 1; j <= numat0; ++j) al[l++] = a[i][j];
  double det = 0.0;
  osinv(al.data(), numat0, det);
  std::vector<double> q(numat + 1, 0.0);
  l = 0;
  for (int i = 1; i <= numat; ++i) {
    double qi = 0.0;
    for (int j = 1; j <= numat0; ++j) qi += al[(i - 1) * numat0 + (j - 1)] * b[j];
    q[i] = qi;
  }
  (void)det;  // dummy use of det
  // AVERAGE CHARGES EQUIVALENT BY SYMMETRY
  if (keywrd.find("SYMAVG") != std::string::npos) {
    std::vector<std::vector<bool>> cequiv(numat + 1, std::vector<bool>(numat + 1, false));
    std::vector<double> qsc(numat + 1, 0.0), qq2(numat + 1, 0.0);
    for (int i = 1; i <= numat; ++i)
      for (int j = 1; j <= numat; ++j) {
        cequiv[i][j] = false;
        if (std::fabs(std::fabs(ch[i]) - std::fabs(ch[j])) >= 1e-5) continue;
        cequiv[i][j] = true;
      }
    for (int i = 1; i <= numat; ++i) {
      int ieq = 0;
      qsc[i] = 0.0;
      for (int j = 1; j <= numat; ++j) {
        if (!cequiv[i][j]) continue;
        qsc[i] = qsc[i] + std::fabs(q[j]);
        ieq = ieq + 1;
      }
      qq2[i] = q[i] / std::fabs(q[i]) * qsc[i] / ieq;
    }
    for (int i = 1; i <= numat; ++i) q[i] = qq2[i];
  }
  // ROOT MEAN SQUARE (RMS) AND RELATIVE RMS
  rms = 0.0;
  rrms = 0.0;
  for (int i = 1; i <= nmep; ++i) {
    double epi = ep[i];
    double epc = 0.0;
    for (int k = 1; k <= numat; ++k) epc += q[k] * d1[k][i];
    epc = epc - epi;
    rms = rms + epc * epc;
    rrms = rrms + epi * epi;
  }
  rms = std::sqrt(rms / nmep);
  rrms = rms / std::sqrt(rrms / nmep);
  rms = rms * au;

  std::printf("    ATOM NO.    TYPE    CHARGE\n");
  for (int i = 1; i <= numat; ++i)
    std::printf("%5d%9s%2s%10.4f\n", i, "", elemnt[nat[i]].c_str(), q[i]);

  std::printf("\n  NUMBER OF POINTS %6d\n", nmep);
  std::printf("  RMS DEVIATION %9.4f\n", rms);
  std::printf("  RRMS DEVIATION %9.4f\n", rrms);

  if (iz != 0) return;
  // DIPOLE MOMENT FOR NEUTRAL MOLECULE
  {
    std::vector<double> rx(numat, 0.0), ry(numat, 0.0), rz(numat, 0.0);
    for (int i = 0; i < numat; ++i) {
      rx[i] = co[1][i + 1];
      ry[i] = co[2][i + 1];
      rz[i] = co[3][i + 1];
    }
    dipx = dipx + ddot(numat, rx, 1, q, 1);
    dipy = dipy + ddot(numat, ry, 1, q, 1);
    dipz = dipz + ddot(numat, rz, 1, q, 1);
  }
  dipx = dipx * bohr1;
  dipy = dipy * bohr1;
  dipz = dipz * bohr1;
  double dip = std::sqrt(dipx * dipx + dipy * dipy + dipz * dipz);
  std::printf("\n\n DIPOLE MOMENT EVALUATED FROM THE MEP CHARGES\n");
  std::printf("       D(X)     D(Y)     D(Z)     TOTAL\n");
  std::printf("%3s%9.4f%9.4f%9.4f%9.4f\n", "", dipx * cf, dipy * cf, dipz * cf, dip * cf);
}

void mepmap(std::vector<std::vector<double>>& c1) {
  std::vector<std::vector<double>> t(4, std::vector<double>(4, 0.0));
  t[1][1] = t[2][2] = t[3][3] = 1.0;
  double step = 0.1;
  static int need1 = 0, natm = 0, icase = 0, iback = 0;
  std::vector<double> r0(4, 0.0);
  double cut = 999.0, z0 = 0.0;
  std::printf("\n\n\n           MOLECULAR ELECTROSTATIC POTENTIAL\n\n"
              " REFERENCES  B.WANG AND G.P.FORD J.COMP.CHEM. 15, 200:207 (1994)\n"
              "            G.P.FORD AND B.WANG J.COMPT.CHEM. 14(1993)1101\n");
  double xm = 0.0, ym = 0.0;
  meprot(coord, r0, t, xm, ym, icase, step, c1, z0, iback);
  if (moperr) return;
  if (iback > 0) return;
  bool prtmep = keywrd.find(" PRTMEP") != std::string::npos;
  bool minmep = keywrd.find(" MINMEP") != std::string::npos;
  std::vector<std::vector<double>> cw(4, std::vector<double>(numat + 1, cut));
  std::vector<double> vmin(numat + 1, cut);
  std::vector<int> mx(numat + 1, 0), my(numat + 1, 0);

  int nonzo = 0;
  std::vector<double> pp(static_cast<std::size_t>(lm61) + 1, 0.0);
  packp(p, pp, nonzo);

  int npx = static_cast<int>(xm / step) * 2 + 1;
  int npy;
  if (icase == 2 || icase == 3) {
    ym = ym + 3.0;
    npy = static_cast<int>(ym / step) + 1;
  } else {
    npy = static_cast<int>(ym / step) * 2 + 1;
  }
  FILE* fp = nullptr;
  if (prtmep) {
    fp = fopen(mep_fn.empty() ? "mep.dat" : mep_fn.c_str(), "a");
    if (fp) {
      std::fprintf(fp, "%3d%5d%5d%8.2f%8.2f%8.2f\n\n", numat, npx, npy, -xm, -ym, step);
      std::fprintf(fp, "%22s\n", "INSTRUCTION");
      std::fprintf(fp, "1) NX*NY = %3d*%3d =%6d GRIDS;\n", npx, npy, npx * npy);
      std::fprintf(fp, "2) NY VALUES IN EACH OF THE NX BLOCKS;\n");
      std::fprintf(fp, "3) X,Y COORDINATES (Z=%6.2f ) OF THE GRIDS:\n", z0);
      std::fprintf(fp, "   X(I)=%8.2f+%6.2f*(I-1)    I=1,2,...,NX\n", -xm, step);
      std::fprintf(fp, "   Y(J)=%8.2f+%6.2f*(J-1)    I=1,2,...,NY\n", -ym, step);
      std::fprintf(fp, "4) THE MOLECULAR COORDINATES (X,Y) ON THE 2D-MEP CONTOUR MAP\n");
      for (int i = 1; i <= numat; ++i)
        std::fprintf(fp, "%9d%4d%12.5f%12.5f\n", i, nat[i], c1[1][i], c1[2][i]);
      std::fprintf(fp, "\n");
    }
  }
  // CREATE THE POINTS AND CALCULATE MEP
  int ip = 0;
  double t11 = t[1][1], t12 = t[1][2], t13 = t[1][3];
  double t21 = t[2][1], t22 = t[2][2], t23 = t[2][3];
  double t31 = t[3][1], t32 = t[3][2], t33 = t[3][3];
  double xx0 = z0 * t31, yy0 = z0 * t32, zz0 = z0 * t33;
  double xxw = xx0 + r0[1], yyw = yy0 + r0[2], zzw = zz0 + r0[3];
  double step1 = step;
  std::vector<double> uu(201, cut);
  std::vector<double> w(4, 0.0), ria(numat + 1, 0.0);
  int imin = 1;
l40:
  double xx1 = xx0 + r0[1], yy1 = yy0 + r0[2], zz1 = zz0 + r0[3];
  for (int m = 1; m <= npx; ++m) {
    double stepm = (m - 1) * step1;
    double xx2 = xx1 + stepm * t11;
    double yy2 = yy1 + stepm * t12;
    double zz2 = zz1 + stepm * t13;
    for (int n = 1; n <= npy; ++n) {
      uu[n] = cut;
      double stepn = (n - 1) * step1;
      w[1] = xx2 + stepn * t21;
      w[2] = yy2 + stepn * t22;
      w[3] = zz2 + stepn * t23;
      if (need1 == 0) {
        // CHECK IF TOO CLOSE TO NUCLEI
        double rmin = 81.0;
        bool skip = false;
        for (int i = 1; i <= numat; ++i) {
          double rr = (w[1] - c1[1][i]) * (w[1] - c1[1][i]) +
                      (w[2] - c1[2][i]) * (w[2] - c1[2][i]) +
                      (w[3] - c1[3][i]) * (w[3] - c1[3][i]);
          if (rr < 0.25) { skip = true; break; }
          if (nat[i] > 1 && rr < 0.36) { skip = true; break; }
          ria[i] = std::sqrt(rr);
          if (rmin <= rr) continue;
          rmin = rr;
          imin = i;
        }
        if (skip) continue;
        double ui = 0.0;
        pmepco(pp, ria, w, ui, c1, nonzo, 1);
        uu[n] = ui;
        if (ui < vmin[imin]) {
          vmin[imin] = ui;
          cw[1][imin] = w[1]; cw[2][imin] = w[2]; cw[3][imin] = w[3];
          mx[imin] = m;
          my[imin] = n;
        }
      } else {
        double ui = 0.0;
        pmepco(pp, ria, w, ui, c1, nonzo, 0);
        if (ui < vmin[natm]) {
          vmin[natm] = ui;
          cw[1][natm] = w[1]; cw[2][natm] = w[2]; cw[3][natm] = w[3];
        }
      }
      ip = ip + 1;
    }
    if (!(need1 == 0 && prtmep)) continue;
    if (fp) {
      for (int ny = 1; ny <= npy; ++ny) std::fprintf(fp, "%11.5f", uu[ny]);
      std::fprintf(fp, "\n");
    }
  }

  if (need1 == 0) std::printf("\n GRIDS NX*NY=%3d *%3d\n", npx, npy);
  if (prtmep && need1 == 0) std::printf("=> CONTOUR FILE IS PRINTED ON CHANNEL 7\n");
  if (need1 <= 0) {
    // REMOVE THE FALSE MINIMA
    int minima = 0;
    for (int i = 1; i <= numat; ++i) {
      if (vmin[i] > 0.0) continue;
      bool keep = true;
      for (int m = mx[i] - 1; m <= mx[i] + 1; ++m) {
        double stepm = (m - 1) * step;
        double xx2 = xxw + stepm * t11;
        double yy2 = yyw + stepm * t12;
        double zz2 = zzw + stepm * t13;
        for (int n = my[i] - 1; n <= my[i] + 1; ++n) {
          double stepn = (n - 1) * step;
          w[1] = xx2 + stepn * t21;
          w[2] = yy2 + stepn * t22;
          w[3] = zz2 + stepn * t23;
          ip = ip + 1;
          double ui = 0.0;
          pmepco(pp, ria, w, ui, c1, nonzo, 0);
          if (ui >= vmin[i] - 0.001) continue;
          vmin[i] = 99.0;
          keep = false;
          break;
        }
        if (!keep) break;
      }
      if (!keep) continue;
      minima = minima + 1;
    }
    if (minima == 0) {
      std::printf(" NO MEP MINIMUM FOUND IN THIS PLANE\n");
      if (fp) fclose(fp);
      return;
    }
    if (!minmep) goto l110;
    need1 = 1;
    step1 = step * 0.2;
    npx = 11;
    npy = 11;
  }
  natm = natm + 1;
  if (natm > numat) goto l110;
  while (vmin[natm] > 0.0) {
    natm = natm + 1;
    if (natm > numat) goto l110;
  }
  r0[1] = cw[1][natm] - step * (t11 + t21);
  r0[2] = cw[2][natm] - step * (t12 + t22);
  r0[3] = cw[3][natm] - step * (t13 + t23);
  goto l40;
l110:
  if (need1 > 0) std::printf("TOTAL NUMBER OF THE POINTS CALCULATED %8d\n", ip);
  std::printf("\n MINIMA IN THE SPECIFIED PLANE\n"
              "     ATOM       UMIN      DISTANCE          XYZ COORDINATE\n"
              "             KCAL/MOL     ANG           IN MOL COORD\n");
  for (int i = 1; i <= numat; ++i) {
    if (vmin[i] >= 0.0) continue;
    double dd = std::sqrt((c1[1][i] - cw[1][i]) * (c1[1][i] - cw[1][i]) +
                          (c1[2][i] - cw[2][i]) * (c1[2][i] - cw[2][i]) +
                          (c1[3][i] - cw[3][i]) * (c1[3][i] - cw[3][i]));
    std::printf("%3d%2s%2s%9.1f%4s%6.2f%4s%9.3f%9.3f%9.3f\n", i, "", elemnt[nat[i]].c_str(),
                vmin[i], "", dd, "", cw[1][i], cw[2][i], cw[3][i]);
  }
  if (fp) fclose(fp);
}

void meprot(const std::vector<std::vector<double>>& c, std::vector<double>& r0,
            std::vector<std::vector<double>>& t, double& xm, double& ym,
            int& icase, double step, std::vector<std::vector<double>>& c1,
            double& z0, int& iback) {
  int n1 = 1, n2 = 1, n3 = 1;
  if (keywrd.find(" PMEPR") != std::string::npos) {
    // F90 reads a data card from unit ir; C++ reads one line from stdin.
    if (!std::getline(std::cin, line)) {
      iback = 1;
      return;
    }
    line = upcase(line);
    std::size_t i = line.find("ICASE");
    if (i != std::string::npos) {
      icase = static_cast<int>(std::lround(reada(line, static_cast<int>(i))));
      i = line.find(" N1");
      n1 = static_cast<int>(std::lround(reada(line, static_cast<int>(i + 4))));
      i = line.find(" N2");
      n2 = static_cast<int>(std::lround(reada(line, static_cast<int>(i + 4))));
      i = line.find(" N3");
      n3 = static_cast<int>(std::lround(reada(line, static_cast<int>(i + 4))));
      i = line.find(" Z0");
      z0 = reada(line, static_cast<int>(i + 4));
    } else {
      std::sscanf(line.c_str(), "%d %d %d %d %lf", &icase, &n1, &n2, &n3, &z0);
    }
    std::printf("\n INPUT CARD FOR PMEPR: ICASE,I,J,K,Z0 =%4d%4d%4d%4d%8.2f\n",
                icase, n1, n2, n3, z0);
    if (icase < 1 || icase > 3) {
      std::printf(" ICASE SHOULD BE 1,2, OR 3 AS PMEPR IS SPECIFIED\n PMEPR WAS NOT EXECUTED\n");
      iback = 1;
      return;
    }
    if (n1 == n2 || n1 == n3 || n2 == n3 || n1 > numat || n2 > numat || n3 > numat) {
      std::printf(" THE THREE REFERENTIAL ATOMS ARE UNREASONABLE, \n PMEPR WAS NOT EXECUTED\n");
      iback = 1;
      return;
    }
    if (icase == 1) {
      if (std::fabs(z0) < 0.0001)
        std::printf("MEPS ARE IN THE PLANE DEFINED BY ATOMS%3d%3d%3d\n", n1, n2, n3);
      else
        std::printf("MEPS ARE IN THE PLANE%6.2f ABOVE THAT DEFINED BY ATOMS%3d%3d%3d\n", z0, n1, n2, n3);
    }
    if (icase == 2) {
      if (std::fabs(z0) < 0.0001)
        std::printf("MEPS ARE IN THE PLANE THROUGH ATOMS%3d%3d%3d PERPENDICULAR TO THE LINE %3d-%3d\n",
                    n1, n2, n3, n1, n2);
      else
        std::printf("MEPS ARE IN THE PLANE THROUGH ATOMS%3d%3d%3d PERPENDICULAR TO %3d-%3d AT Z0=%6.2f\n",
                    n1, n2, n3, n1, n2, z0);
    }
    if (icase == 3) {
      if (std::fabs(z0) < 0.0001)
        std::printf("MEPS ARE IN THE PLANE BISECTING THE ANGLE (%3d%3d%3d)\n", n2, n1, n3);
      else
        std::printf("MEPS ARE IN THE PLANE BISECTING THE ANGLE (%3d%3d%3d)\n"
                    "     BUT NON-ZERO VALUE Z0=%6.2f MAKES NO SENSE AND IS IGNORED\n",
                    n2, n1, n3, z0);
    }
  }
  if (icase <= 1) {
    for (int i = 1; i <= numat; ++i)
      for (int d = 1; d <= 3; ++d) c1[d][i] = c[d][i];
    goto l150;
  }
  // MOVE THE ORIGIN TO ATOM N1
  for (int i = 1; i <= numat; ++i)
    for (int d = 1; d <= 3; ++d) c1[d][i] = c[d][i] - c[d][n1];

  double x2 = c1[1][n2], y2 = c1[2][n2], z2 = c1[3][n2];
  double rr2 = x2 * x2 + y2 * y2 + z2 * z2;
  double r2 = std::sqrt(rr2);
  double r21 = 1.0 / r2;
  double x3 = c1[1][n3], y3 = c1[2][n3], z3 = c1[3][n3];
  double rr3 = x3 * x3 + y3 * y3 + z3 * z3;
  double r3 = std::sqrt(rr3);
  double r31 = 1.0 / r3;
  double cos1 = (rr3 + rr2 - (x3 - x2) * (x3 - x2) - (y3 - y2) * (y3 - y2) -
                 (z3 - z2) * (z3 - z2)) * (0.5 * r31 * r21);
  if (1.0 - std::fabs(cos1) < 0.01) {
    std::printf("N2-N1-N3=%9.3f DEGREES, ALMOST LINEAR!\n", std::acos(cos1));
    mopend("Error in PMEP");
    return;
  }

  if (icase == 1 || icase == 2) {
    // CHOOSE X-AXIS FOR CASE 1,2
    t[1][1] = x2 * r21;
    t[1][2] = y2 * r21;
    t[1][3] = z2 * r21;
    // Y-AXIS FOR CASE 1, Z CASE 2
    double p41 = r3 * cos1 * r21;
    double x34 = x3 - p41 * x2;
    double y34 = y3 - p41 * y2;
    double z34 = z3 - p41 * z2;
    double r341 = 1.0 / std::sqrt(x34 * x34 + y34 * y34 + z34 * z34);
    if (icase == 1) {
      t[2][1] = x34 * r341;
      t[2][2] = y34 * r341;
      t[2][3] = z34 * r341;
    } else {
      t[3][1] = -x34 * r341;
      t[3][2] = -y34 * r341;
      t[3][3] = -z34 * r341;
    }
  }
  // CHOOSE X, Z-AXES FOR CASE 3
  if (icase == 3) {
    double x4, y4, z4;
    if (r3 > r2) {
      double x6 = r2 * r31;
      double y6 = x6 * y3;
      double z6 = x6 * z3;
      x6 = x6 * x3;
      double x62 = x6 - x2, y62 = y6 - y2, z62 = z6 - z2;
      double r621 = 1.0 / std::sqrt(x62 * x62 + y62 * y62 + z62 * z62);
      t[3][1] = -x62 * r621;
      t[3][2] = -y62 * r621;
      t[3][3] = -z62 * r621;
      x4 = (x2 + x6) * 0.5;
      y4 = (y2 + y6) * 0.5;
      z4 = (z2 + z6) * 0.5;
    } else {
      double x6 = r3 * r21;
      double y6 = x6 * y2;
      double z6 = x6 * z3;
      x6 = x6 * x2;
      double x63 = x6 - x3, y63 = y6 - y3, z63 = z6 - z3;
      double r631 = 1.0 / std::sqrt(x63 * x63 + y63 * y63 + z63 * z63);
      t[3][1] = x63 * r631;
      t[3][2] = y63 * r631;
      t[3][3] = z63 * r631;
      x4 = (x3 + x6) * 0.5;
      y4 = (y3 + y6) * 0.5;
      z4 = (z3 + z6) * 0.5;
    }
    double r41 = 1.0 / std::sqrt(x4 * x4 + y4 * y4 + z4 * z4);
    t[1][1] = x4 * r41;
    t[1][2] = y4 * r41;
    t[1][3] = z4 * r41;
  }
  // Z FOR ICASE=1, Y FOR ICASE=2 & 3
  double x5 = 0.0, y5 = 0.0;
  double z5 = 1.0;
  if (std::fabs(x2) > 0.1) {
    y5 = -(x2 * z3 - x3 * z2) / (x2 * y3 - x3 * y2);
    x5 = -(y5 * y2 + z2) / x2;
    goto l130;
  }
  if (std::fabs(x3) > 0.1) {
    y5 = -(x3 * z2 - x2 * z3) / (x3 * y2 - x2 * y3);
    x5 = -(y3 * y5 + z3) / x3;
    goto l130;
  }
  if (std::fabs(y3) > 0.1) {
    x5 = -(y3 * z2 - y2 * z3) / (x2 * y3 - x3 * y2);
    y5 = -(x5 * x3 + z3) / y3;
    goto l130;
  }
  if (std::fabs(y2) > 0.1) {
    x5 = -(y2 * z3 - y3 * z2) / (x3 * y2 - x2 * y3);
    y5 = -(x5 * x2 + z2) / y2;
    goto l130;
  } else {
    std::printf("X21,Y21,X31,Y31=%8.4f%8.4f%8.4f%8.4f WHY THEY ARE SO SMALL?\n", x2, y2, x3, y3);
    mopend("Error in PMEP");
    return;
  }
l130: {
  double r51 = 1.0 / std::sqrt(x5 * x5 + y5 * y5 + 1.0);
  int k = 3;
  if (icase != 1) k = 2;
  t[k][1] = x5 * r51;
  t[k][2] = y5 * r51;
  t[k][3] = z5 * r51;
  // ROTATE THE MOLECULE TO NEW XYZ COORDINATES
  for (int kk = 1; kk <= numat; ++kk) {
    double xi = t[1][1] * c1[1][kk] + t[1][2] * c1[2][kk] + t[1][3] * c1[3][kk];
    double yi = t[2][1] * c1[1][kk] + t[2][2] * c1[2][kk] + t[2][3] * c1[3][kk];
    c1[3][kk] = t[3][1] * c1[1][kk] + t[3][2] * c1[2][kk] + t[3][3] * c1[3][kk];
    c1[1][kk] = xi;
    c1[2][kk] = yi;
  }
}
  // NON-WEIGHTED CENTER
l150: {
  double c0[4] = {0, 0, 0, 0};
  for (int j = 1; j <= 2; ++j) {
    double cs = 1e6, cl = -1e6;
    int is = 1, il = 1;
    for (int i = 1; i <= numat; ++i) {
      if (c1[j][i] < cs) { cs = c1[j][i]; is = i; }
      if (c1[j][i] <= cl) continue;
      cl = c1[j][i];
      il = i;
    }
    c0[j] = (cs + cl) / 2.0;
    r0[j] = cl - c0[j] + 2.0;
    if (nat[is] <= 1 && nat[il] <= 1) continue;
    r0[j] = r0[j] + 0.4;
  }
  // MOVE TO THE NEW CENTER
  for (int i = 1; i <= numat; ++i) {
    c1[1][i] = c1[1][i] - c0[1];
    c1[2][i] = c1[2][i] - c0[2];
  }
  xm = static_cast<int>(r0[1] / step) * step;
  ym = static_cast<int>(r0[2] / step) * step;
  if (icase < 1 || icase > 3) {
    r0[1] = ((-xm) + c0[1]) * t[1][1] + ((-ym) + c0[2]) * t[2][1];
    r0[2] = ((-xm) + c0[1]) * t[1][2] + ((-ym) + c0[2]) * t[2][2];
    r0[3] = ((-xm) + c0[1]) * t[1][3] + ((-ym) + c0[2]) * t[2][3];
  } else {
    if (icase == 1) {
      r0[1] = ((-xm) + c0[1]) * t[1][1] + ((-ym) + c0[2]) * t[2][1] + c[1][n1];
      r0[2] = ((-xm) + c0[1]) * t[1][2] + ((-ym) + c0[2]) * t[2][2] + c[2][n1];
      r0[3] = ((-xm) + c0[1]) * t[1][3] + ((-ym) + c0[2]) * t[2][3] + c[3][n1];
    } else {
      r0[1] = ((-xm) + c0[1]) * t[1][1] + c0[2] * t[2][1] + c[1][n1];
      r0[2] = ((-xm) + c0[1]) * t[1][2] + c0[2] * t[2][2] + c[2][n1];
      r0[3] = ((-xm) + c0[1]) * t[1][3] + c0[2] * t[2][3] + c[3][n1];
    }
  }
}
}

void packp(const std::vector<double>& p, std::vector<double>& pp, int& mn) {
  mn = 0;
  for (int i = 1; i <= numat; ++i) {
    int ia = nfirst[i];
    int ib = nlast[i];
    for (int j = ia; j <= ib; ++j) {
      int j1st = j * (j - 1) / 2;
      if (j - ia + 1 > 0) {
        for (int t = 0; t < j - ia + 1; ++t) pp[mn + 1 + t] = p[j1st + ia + t];
        mn = j - ia + 1 + mn;
      }
    }
  }
}
