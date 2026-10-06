// static_polarizability.cpp — C++ translation of "static_polarizability.F90"
// (static_polarizability + ffhpol + dipind, 548-line Fortran).
//
// Finite-field calculation of molecular electric response properties:
// dipole moment, polarizability, first and second hyperpolarizability.
// Field steps drive the compfg SCF driver (skeleton) 36 times; results are
// summarized as H.o.F. (heat of formation) and Dipole expansions, then the
// polarizability tensors are diagonalized by the packed Jacobi eigensolver rsp.

#include "static_polarizability.h"
#include "dipind.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "axis.h"
#include "chanel_C.h"
#include "chrge.h"
#include "common_arrays_C.h"
#include "compfg.h"
#include "elemts_C.h"
#include "funcon_C.h"
#include "gmetry.h"
#include "ijbo.h"
#include "matout.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "pol_vol.h"
#include "rsp.h"
#include "vecprt.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace elemts_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

namespace {

// Flatten the 1-based (3 x n) internal coordinate array geo[d][i] into the
// 0-based column-major vector compfg expects: xparam[(i-1)*3+(d-1)] = geo[d][i].
std::vector<double> flatten_geo(const std::vector<std::vector<double>>& g, int n) {
  std::vector<double> f(3 * n, 0.0);
  for (int i = 1; i <= n; ++i)
    for (int d = 1; d <= 3; ++d) f[(i - 1) * 3 + (d - 1)] = g[d][i];
  return f;
}

}  // namespace

void static_polarizability() {
  bool let = (keywrd.find(" LET") != std::string::npos);
  std::printf(" ******************** FINITE-FIELD POLARIZABILITIES ********************\n\n"
              "    THE F-F METHOD IS PERFORMED USING BOTH AN ENERGY\n"
              "    AND DIPOLE MOMENT EXPANSION.  THESE RESULTS ARE\n"
              "    LISTED BELOW AS \"H.o.F.\" AND \"Dipole\", RESPECTIVELY.\n");
  gmetry(geo, coord);

  // Orient the molecule with the moments of inertia.  This is done to ensure
  // a unique, reproduceable set of directions.  If LET is specified, the
  // input orientation will be used.
  if (!let) {
    double sum = 0.0, sum2 = 0.0, sum3 = 0.0;
    double rotvec[4][4] = {};
    axis(sum, sum2, sum3, rotvec);
    std::printf("\n ROTATION MATRIX FOR ORIENTATION OF MOLECULE:\n");
    for (int i = 1; i <= 3; ++i) {
      std::printf("     %12.6f%12.6f%12.6f\n", rotvec[i][1], rotvec[i][2], rotvec[i][3]);
    }
    // ROTATE ATOMS
    for (int i = 1; i <= numat; ++i) {
      for (int j = 1; j <= 3; ++j) {
        double s = 0.0;
        for (int k = 1; k <= 3; ++k) s += coord[k-1][i] * rotvec[k][j];
        geo[j][i] = s;
      }
    }
    for (int i = 1; i <= numat; ++i)
      for (int j = 1; j <= 3; ++j) coord[j-1][i] = geo[j][i];
    std::printf("\n\n          CARTESIAN COORDINATES \n\n");
    std::printf("    NO.       ATOM        X            Y            Z\n\n");
    int l = 0;
    for (int i = 1; i <= numat; ++i) {
      if (nat[i] == 99 || nat[i] == 107) continue;
      ++l;
      std::printf("%6d%8s%4s%10.4f%10.4f%10.4f\n", l, "", elemnt[nat[i]].c_str(),
                  coord[0][l], coord[1][l], coord[2][l]);
    }
  }

  // SET UP THE VARIABLES IN XPARAM AND LOC, THESE ARE IN CARTESIAN COORDINATES.
  numat = 0;
  double sumx = 0.0, sumy = 0.0, sumz = 0.0;
  for (int i = 1; i <= natoms; ++i) {
    if (labels[i] == 99 || labels[i] == 107) continue;
    numat = numat + 1;
    labels[numat] = labels[i];
    sumx = sumx + coord[0][numat];
    sumy = sumy + coord[1][numat];
    sumz = sumz + coord[2][numat];
    for (int d = 1; d <= 3; ++d) geo[d][numat] = coord[d-1][numat];
  }
  sumx = sumx / numat;
  sumy = sumy / numat;
  sumz = sumz / numat;
  double summax = 0.0;
  for (int i = 1; i <= numat; ++i) {
    geo[1][i] = geo[1][i] - sumx;
    summax = std::max(std::fabs(geo[1][i]), summax);
    geo[2][i] = geo[2][i] - sumy;
    summax = std::max(std::fabs(geo[2][i]), summax);
    geo[3][i] = geo[3][i] - sumz;
    summax = std::max(std::fabs(geo[3][i]), summax);
  }

  ndep = 0;
  natoms = numat;
  nvar = 0;
  std::fill(na.begin(), na.end(), 0);
  {
    std::vector<double> grad(1, 0.0);
    compfg(flatten_geo(geo, numat), true, escf, true, grad, false);
  }
  std::printf("\n\n ENERGY OF \"REORIENTED\" SYSTEM WITHOUT FIELD:%15.5f Kcal/mol\n", escf);

  ffhpol();
}

void ffhpol() {
  // Energy: a.u. to kcal/mole
  double autokc = fpc_9 * ev;
  // Dipole: a.u. to debye
  double autodb = a0 * fpc_8 * fpc_1 * 1.0e-10;
  // Electric Field: a.u. to volt/meter
  double autovm = ev / a0;  // = 51.42
  int nbdip = 1, nbcnt = 4, ngcnt = 4, counter = 0;
  double heat0 = escf;
  bool debug = (keywrd.find("DEBUG") != std::string::npos);

  std::vector<double> dipe4(4, 0.0), dipdp(4, 0.0), eigs(4, 0.0), vectrs(10, 0.0);
  double heat1p = 0.0, heat1m = 0.0, heat2p = 0.0, heat2m = 0.0;
  double hpp = 0.0, hpm = 0.0, hmm = 0.0, hmp = 0.0;
  double h2pp = 0.0, h2pm = 0.0, h2mm = 0.0, h2mp = 0.0;
  std::vector<double> apole4(7, 0.0), apoldp(7, 0.0);
  std::vector<double> dip1p(4, 0.0), dip1m(4, 0.0), dip2p(4, 0.0), dip2m(4, 0.0);
  std::vector<double> grad(1, 0.0);
  const char* ch_xyz[4] = {"", "X", "Y", "Z"};

  // FIELD STRENGTH IN A.U.
  double efval = 0.002;
  std::printf("\n\n APPLIED ELECTRIC FIELD MAGNITUDE: %15.5f VOLTS PER ANGSTROM\n", efval * ev / a0);
  double sfe = 1.0 / efval;

  // CALCULATE THE POLARIZABILITY AND HYPERPOLARIZABILITIES ALONG THE THREE
  // PRINCIPLE AXES.
  for (int id = 1; id <= 3; ++id) {
    if (debug) std::printf("\n\n ****** %c DIRECTION *****\n", ch_xyz[id][0]);

    // ZERO THE FIELD
    std::fill(efield, efield + 4, 0.0);
    double hnuc = 0.0;
    for (int i = 1; i <= numat; ++i)
      hnuc += efval * geo[id][i] * tore[nat[i]] * autovm;
    hnuc = hnuc * fpc_9;

    // +E(ID)
    efield[id] = efval;
    compfg(flatten_geo(geo, numat), true, heat1p, true, grad, false);
    ++counter;
    std::printf("    Step%3d of 36 done (+%c)\n", counter, ch_xyz[id][0]);
    dipind(dip1p);
    // -E(ID)
    efield[id] = -efval;
    compfg(flatten_geo(geo, numat), true, heat1m, true, grad, false);
    ++counter;
    std::printf("    Step%3d of 36 done (-%c)\n", counter, ch_xyz[id][0]);
    dipind(dip1m);
    // +2E(ID)
    efield[id] = 2.0 * efval;
    compfg(flatten_geo(geo, numat), true, heat2p, true, grad, false);
    ++counter;
    std::printf("    Step%3d of 36 done (+2%c)\n", counter, ch_xyz[id][0]);
    dipind(dip2p);
    // -2E(ID)
    efield[id] = -2.0 * efval;
    compfg(flatten_geo(geo, numat), true, heat2m, true, grad, false);
    ++counter;
    std::printf("    Step%3d of 36 done (-2%c)\n", counter, ch_xyz[id][0]);
    dipind(dip2m);
    if (debug) {
      std::printf(" ENERGIES AT:            F                    2F\n");
      std::printf("   + %20.10f   %20.10f\n", heat1p, heat2p);
      std::printf("   - %20.10f   %20.10f\n", heat1m, heat2m);
    }

    // DIPOLE
    double eterm = (1.0 / 12.0) * (heat2p - heat2m) - (2.0 / 3.0) * (heat1p - heat1m);
    dipe4[id] = eterm * sfe / autokc;

    // ALPHA
    int ivl = (id * (id + 1)) / 2;
    eterm = 2.5 * heat0 - (4.0 / 3.0) * (heat1p + heat1m) + (1.0 / 12.0) * (heat2p + heat2m);
    apole4[ivl] = eterm * sfe * sfe / autokc;

    // DIPOLE CALCULATIONS
    double dmu = (2.0 / 3.0) * (dip1p[id] + dip1m[id]) - (1.0 / 6.0) * (dip2p[id] + dip2m[id]);
    dipdp[id] = dmu / autodb;
    double ae = (2.0 / 3.0) * (dip1p[id] - dip1m[id]) - (1.0 / 12.0) * (dip2p[id] - dip2m[id]);
    apoldp[ivl] = ae * sfe / autodb;
    for (int kd = 1; kd <= 3; ++kd) {
      if (kd < id) {
        int kvl = (id * (id - 1)) / 2 + kd;
        double aki = (2.0 / 3.0) * (dip1p[kd] - dip1m[kd]) - (1.0 / 12.0) * (dip2p[kd] - dip2m[kd]);
        apoldp[kvl] = aki * sfe / autodb;
      }
      if (kd == id) continue;
      nbdip = nbdip + 1;
    }

    // NOW CALCULATE THE OFF AXIS RESULTS.
    int idm1 = id - 1;
    for (int jd = 1; jd <= idm1; ++jd) {
      double hnucj = 0.0;
      for (int i = 1; i <= numat; ++i)
        hnucj += efval * geo[jd][i] * tore[nat[i]] * ev / a0;
      hnucj = hnucj * fpc_9;
      std::fill(efield, efield + 4, 0.0);

      // DIAGONAL FIELDS WITH COMPONENTS EQUAL TO EFVAL
      efield[id] = efval;
      efield[jd] = efval;
      compfg(flatten_geo(geo, numat), true, hpp, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (+%c+%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      dipind(dip1p);
      efield[jd] = -efval;
      compfg(flatten_geo(geo, numat), true, hpm, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (+%c-%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      dipind(dip1p);
      efield[id] = -efval;
      compfg(flatten_geo(geo, numat), true, hmm, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (-%c-%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      dipind(dip1p);
      efield[jd] = efval;
      compfg(flatten_geo(geo, numat), true, hmp, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (-%c+%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      dipind(dip1p);
      hpp = hpp + hnuc + hnucj;
      hpm = hpm + hnuc - hnucj;
      hmm = hmm - hnuc - hnucj;
      hmp = hmp - hnuc + hnucj;
      if (debug) {
        std::printf("\n             +,+             +,-             -,+             -,-\n");
        std::printf("  E %18.6f%18.6f%18.6f%18.6f\n", hpp, hpm, hmp, hmm);
      }

      // DIAGONAL FIELDS WITH COMPONENTS EQUAL TO 2*EFVAL
      efield[id] = efval * 2.0;
      efield[jd] = efval * 2.0;
      compfg(flatten_geo(geo, numat), true, h2pp, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (+2%c+2%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      efield[jd] = -efval * 2.0;
      compfg(flatten_geo(geo, numat), true, h2pm, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (+2%c-2%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      efield[id] = -efval * 2.0;
      compfg(flatten_geo(geo, numat), true, h2mm, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (-2%c-2%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      efield[jd] = efval * 2.0;
      compfg(flatten_geo(geo, numat), true, h2mp, true, grad, false);
      ++counter;
      std::printf("    Step%3d of 36 done (-2%c+2%c)\n", counter, ch_xyz[id][0], ch_xyz[jd][0]);
      h2pp = h2pp + 2.0 * (hnuc + hnucj);
      h2pm = h2pm + 2.0 * (hnuc - hnucj);
      h2mm = h2mm - 2.0 * (hnuc + hnucj);
      h2mp = h2mp - 2.0 * (hnuc - hnucj);
      if (debug) {
        std::printf(" 2E %18.6f%18.6f%18.6f%18.6f\n", h2pp, h2pm, h2mp, h2mm);
      }

      double aterm = (1.0 / 48.0) * (h2pp - h2pm - h2mp + h2mm) -
                     (1.0 / 3.0) * (hpp - hpm - hmp + hmm);
      double aij = aterm * sfe * sfe / autokc;
      ivl = (id * (id - 1)) / 2 + jd;
      apole4[ivl] = aij;
      nbcnt = nbcnt + 1;
      nbcnt = nbcnt + 1;
      ngcnt = ngcnt + 1;
    }
  }

  // SUMMARIZE THE RESULTS
  std::printf("\n\n %30s DIPOLE %30s\n\n", "******************************", "******************************");
  double dipe4t = std::sqrt(dipe4[1] * dipe4[1] + dipe4[2] * dipe4[2] + dipe4[3] * dipe4[3]);
  double dipe4d = dipe4t * autodb;
  double dipdpt = std::sqrt(dipdp[1] * dipdp[1] + dipdp[2] * dipdp[2] + dipdp[3] * dipdp[3]);
  double dipdpd = dipdpt * autodb;
  std::printf("                     H.o.F.         Dipole\n");
  std::printf("     %c%7s%15.6f%15.6f\n", 'X', "", dipe4[1], dipdp[1]);
  std::printf("     %c%7s%15.6f%15.6f\n", 'Y', "", dipe4[2], dipdp[2]);
  std::printf("     %c%7s%15.6f%15.6f\n", 'X', "", dipe4[3], dipdp[3]);
  std::printf("\n\n MAGNITUDE:  %15.6f%15.6f  (A.U.)\n %12s%15.6f%15.6f  (DEBYE)\n",
              dipe4t, dipdpt, "", dipe4d, dipdpd);

  // FIND EIGENVALUES AND EIGENVECTORS OF POLARIZATION MATRIX.
  std::printf("\n\n %30s POLARIZABILITY %20s\n\n", "******************************", "********************");
  std::printf("\n Polarizability Tensor from Heat of Formation:\n");
  apole4[1] = pol_vol(apole4[1]);
  apole4[3] = pol_vol(apole4[3]);
  apole4[6] = pol_vol(apole4[6]);
  int i3 = 3;
  // apole4 is 1-based (index 0 is a pad); vecprt/rsp consume 0-based buffers.
  vecprt(apole4.data() + 1, (-i3));
  rsp(apole4.data() + 1, i3, eigs.data() + 1, vectrs.data());
  {
    int nr = 3;
    matout(vectrs.data(), eigs.data() + 1, i3, nr, i3);
  }
  double avgpe4 = (eigs[1] + eigs[2] + eigs[3]) / 3.0;
  double avga3 = avgpe4 * 0.14818;
  double avgesu = avgpe4 * 0.296352e-24;
  std::printf("\n Polarizability Tensor from Dipole Moment:\n");
  apoldp[1] = pol_vol(apoldp[1]);
  apoldp[3] = pol_vol(apoldp[3]);
  apoldp[6] = pol_vol(apoldp[6]);
  vecprt(apoldp.data() + 1, (-i3));
  rsp(apoldp.data() + 1, i3, eigs.data() + 1, vectrs.data());
  {
    int nr = 3;
    matout(vectrs.data(), eigs.data() + 1, i3, nr, i3);
  }
  double avgpdp = (eigs[1] + eigs[2] + eigs[3]) / 3.0;
  double avga3d = avgpdp * 0.14818;
  double avgesd = avgpdp * 0.296352e-24;
  std::printf("\n\n Average Polarizability from:     H.o.F         Dipole\n"
              " %24s%15.6f%15.6f  A.U.\n"
              " %24s%15.6f%15.6f  ANG.**3\n"
              " %24s%15.6e%15.6e  ESU\n",
              "", avgpe4, avgpdp, "", avga3, avga3d, "", avgesu, avgesd);

  (void)nbdip;
  (void)nbcnt;
  (void)ngcnt;
}
