// react1.cpp — C++ translation of MOPAC 2016 "react1.F90".
//
//  REACT1 DETERMINES THE TRANSITION STATE OF A CHEMICAL REACTION.
//   REACT WORKS BY USING TWO SYSTEMS SIMULTANEOUSLY, THE HEATS OF
//   FORMATION OF BOTH ARE CALCULATED, THEN THE MORE STABLE ONE
//   IS MOVED IN THE DIRECTION OF THE OTHER. AFTER A STEP THE
//   ENERGIES ARE COMPARED, AND THE NOW LOWER-ENERGY FORM IS MOVED
//   IN THE DIRECTION OF THE HIGHER-ENERGY FORM. THIS IS REPEATED
//   UNTIL THE SADDLE POINT IS REACHED.
//
//  DOCK rotates and translates the geometry in GEOA so that the root
//       mean square difference between GEOA and GEO is a minimum.
#include "react1.h"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "reada.h"
#include "getgeo.h"
#include "mopend.h"
#include "second.h"
#include "geout.h"
#include "compfg.h"
#include "writmo.h"
#include "ef.h"
#include "flepo.h"
#include "to_screen.h"
#include "xyzint.h"
#include "gmetry.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "elemts_C.h"
#include "common_arrays_C.h"

namespace {

// BLAS ddot over Fortran 1-based arrays.
double ddot(int n, const std::vector<double>& x, int incx,
            const std::vector<double>& y, int incy) {
  double s = 0.0;
  int xi = 1, yi = 1;
  for (int i = 0; i < n; ++i) {
    s += x[xi] * y[yi];
    xi += incx;
    yi += incy;
  }
  return s;
}

// Fortran index(s, sub): 1-based position or 0.
int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// Fortran trim: strip trailing blanks.
std::string trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? std::string() : s.substr(0, e + 1);
}

// Fortran nint: nearest integer.
long nint(double x) { return std::lround(x); }

}  // namespace

void react1() {
  using namespace molkst_C;
  using namespace common_arrays_C;

  int linear = 0, iflag = 1, i, j, k, iloop, maxstp, maxcyc;
  bool gok[2];
  bool gradnt, intl = false, lef;
  double stepmx, x, sumx, sumy, sumz, sum, step0, one, dell, eold;
  double time1, swap, funct1, time2, gold, c1, c2, dist;

  gok[0] = false;
  gok[1] = false;
  lef = index1(keywrd, " EF") != 0;
  gradnt = index1(keywrd, "GRAD") != 0;
  i = index1(keywrd, " BAR");
  stepmx = 0.01;
  if (i != 0) stepmx = reada(keywrd, i);
  maxstp = 10000;

  std::vector<std::vector<int>> lopt(4, std::vector<int>(natoms + 1, 0));
  std::vector<double> xold(3 * numat + 1, 0.0), grold(3 * numat + 1, 0.0);
  std::vector<double> pastor, pbstor;
  std::vector<double> xyz;  // column-major (3 x numat) buffer for xyzint

  if (index1(keywrd, "GEO_REF") == 0) {
    // READ IN THE SECOND GEOMETRY.
    if ((int)geoa.size() < 4) geoa.assign(4, std::vector<double>(natoms + 1, 0.0));
    else for (auto& r : geoa) if ((int)r.size() < natoms + 1) r.resize(natoms + 1, 0.0);
    i = numat;
    getgeo(chanel_C::ir, labels, geoa, coord, lopt, na, nb, nc, intl);
    if (numat != i) {
      // Fortran: rewind(ir) — no-op in the C++ port (input is std::cin /
      // redirected file; getgeo already consumed the failing geometry).
      std::fprintf(stdout, "%s%4d\n", " Number of atoms in the first geometry: ", i);
      std::fprintf(stdout, "%s%4d\n",
                   " Number of atoms in the second geometry:", numat);
      std::fprintf(stdout, "%s\n%s\n",
                   " Second geometry has different number of atoms from first geometry",
                   " Data set suppled:");
      for (i = 1; i <= 100000; ++i) {
        if (!std::getline(std::cin, line)) break;   // iostat < 0 -> EOF
        std::fprintf(stdout, "%s\n", trim(line).c_str());
      }
      mopend("Geometry supplied is faulty - see output file");
      // Dummy use of intl
      if (intl || lopt[1][1] != 0) return;
      return;
    }
  } else {
    gmetry(geoa, coord);
  }
  maxcyc = 100000;
  if (index1(keywrd, " BIGCYCLES") != 0)
    maxcyc = (int)nint(reada(keywrd, index1(keywrd, " BIGCYCLES")));

  // SWAP FIRST AND SECOND GEOMETRIES AROUND SO THAT GEOUT CAN OUTPUT DATA
  // ON SECOND GEOMETRY.
  xyz.assign(3 * numat, 0.0);
  for (i = 1; i <= numat; ++i)
    for (k = 1; k <= 3; ++k) xyz[(i - 1) * 3 + (k - 1)] = geo[k][i];
  for (i = 1; i <= numat; ++i)
    for (k = 1; k <= 3; ++k) { geo[k][i] = geoa[k][i]; geoa[k][i] = xyz[(i - 1) * 3 + (k - 1)]; }
  sumx = 0.0; sumy = 0.0; sumz = 0.0;
  for (j = 1; j <= numat; ++j) {
    sumx = sumx + coord[0][j];
    sumy = sumy + coord[1][j];
    sumz = sumz + coord[2][j];
  }
  sumx = sumx / numat; sumy = sumy / numat; sumz = sumz / numat;
  for (j = 1; j <= numat; ++j) {
    geo[1][j] = coord[0][j] - sumx;
    geo[2][j] = coord[1][j] - sumy;
    geo[3][j] = coord[2][j] - sumz;
  }
  std::fprintf(stdout, "\n\n  CARTESIAN GEOMETRY OF FIRST SYSTEM\n\n");
  for (i = 1; i <= numat; ++i)
    std::fprintf(stdout, "%4d   %.2s   %14.5f%14.5f%14.5f\n", i,
                 elemts_C::elemnt[nat[i]].c_str(), geo[1][i], geo[2][i], geo[3][i]);
  dock(geoa, geo, dist);
  std::fprintf(stdout, "\n\n  CARTESIAN GEOMETRY OF SECOND SYSTEM\n\n");
  for (i = 1; i <= numat; ++i)
    std::fprintf(stdout, "%4d   %.2s   %14.5f%14.5f%14.5f\n", i,
                 elemts_C::elemnt[nat[i]].c_str(), geoa[1][i], geoa[2][i], geoa[3][i]);
  std::fprintf(stdout, "\n\n   \"DISTANCE\":%13.6f\n", dist);
  std::fprintf(stdout, "\n\n  REACTION COORDINATE VECTOR\n\n");
  for (i = 1; i <= numat; ++i)
    std::fprintf(stdout, "%4d   %.2s   %14.5f%14.5f%14.5f\n", i,
                 elemts_C::elemnt[nat[i]].c_str(), geo[1][i] - geoa[1][i],
                 geo[2][i] - geoa[2][i], geo[3][i] - geoa[3][i]);

  // XPARAM HOLDS THE VARIABLE PARAMETERS FOR GEOMETRY IN GEO
  // XOLD   HOLDS THE VARIABLE PARAMETERS FOR GEOMETRY IN GEOA
  if ((int)xparam.size() < 3 * numat + 1) xparam.resize(3 * numat + 1, 0.0);
  sum = 0.0;
  i = 0;
  for (j = 1; j <= numat; ++j) {
    for (k = 1; k <= 3; ++k) {
      i = i + 1;
      grold[i] = 1.0;
      xparam[i] = geo[k][j];
      xold[i] = geoa[k][j];
      sum = sum + (xparam[i] - xold[i]) * (xparam[i] - xold[i]);
    }
  }
  step0 = std::sqrt(sum);
  if (step0 < 1.0e-2) {
    std::fprintf(stdout, "\n\n\n            %s\n            %s\n            %s\n",
                 " THE TWO GEOMETRIES ARE IDENTICAL OR ALMOST IDENTICAL",
                 " A SADDLE CALCULATION INVOLVES A REACTANT AND A PRODUCT",
                 " THESE MUST BE DIFFERENT GEOMETRIES");
    mopend("THE TWO GEOMETRIES IN SADDLE ARE IDENTICAL.");
    return;
  }
  time0 = second(2);
  maxcyc = 100000;
  if (index1(keywrd, " BIGCYCLES") != 0)
    maxcyc = (int)nint(reada(keywrd, index1(keywrd, " BIGCYCLES")));

  one = 1.0;
  dell = 0.1;
  eold = -2000.0;
  time1 = second(2);
  swap = 0.0;
  if ((int)grad.size() < nvar + 1) grad.resize(nvar + 1, 0.0);
  for (iloop = 1; iloop <= maxstp; ++iloop) {
    if (iloop >= maxcyc) tleft = -100.0;
    std::fprintf(stdout, " ");
    for (k = 0; k < 40; ++k) std::fprintf(stdout, "*+");
    std::fprintf(stdout, "\n");

    // THIS METHOD OF CALCULATING 'STEP' IS QUITE ARBITARY, AND NEEDS
    // TO BE IMPROVED BY INTELLIGENT GUESSWORK!
    gnorm = std::max(1.0e-3, gnorm);
    {
      double m1 = std::min(std::min(swap, 0.5), 6.0 / gnorm);
      double m2 = std::min(dell, stepmx * step0);
      step = std::min(0.2, std::min(m1, m2) / step0) * step0;
    }
    swap = swap + 1.0;
    dell = dell + 0.1;
    std::fprintf(stdout, "            BAR SHORTENED BY%12.3f PERCENT\n",
                 step / step0 * 100.0);
    step0 = step0 - step;
    if (step0 < 0.01) break;
    step = step0;
    if (lef) {
      ef(xparam, escf);
      if (moperr) return;
    } else {
      flepo(xparam, nvar, escf);
      if (moperr) return;
    }
    if (linear == 0 && !pa.empty()) {
      linear = mpack;
      pastor = pa;
      pbstor = pb;
    }
    i = 0;
    for (j = 1; j <= numat; ++j)
      for (k = 1; k <= 3; ++k) {
        i = i + 1;
        xparam[i] = geo[k][j];
      }
    if (iflag == 1)
      std::fprintf(stdout, "\n\n          FOR POINT%5d SECOND STRUCTURE\n", iloop);
    else
      std::fprintf(stdout, "\n\n          FOR POINT%5d FIRST  STRUCTURE\n", iloop);
    std::fprintf(stdout, " DISTANCE A - B  %12.6f\n", step);
    sum = step;

    // NOW TO CALCULATE THE "CORRECT" GRADIENTS, SWITCH OFF 'STEP'.
    step = 0.0;
    for (i = 1; i <= nvar; ++i) grad[i] = grold[i];
    compfg(xparam, true, funct1, true, grad, true);
    if (moperr) return;
    for (i = 1; i <= nvar; ++i) grold[i] = grad[i];
    if (gradnt) {
      std::fprintf(stdout, "  ACTUAL GRADIENTS OF THIS POINT\n");
      for (i = 1; i <= nvar; ++i) {
        std::fprintf(stdout, "%10.4f", grad[i]);
        if (i % 8 == 0) std::fprintf(stdout, "\n");
      }
      if (nvar % 8 != 0) std::fprintf(stdout, "\n");
    }
    std::fprintf(stdout, " HEAT            %12.6f\n", funct1);
    gnorm = std::sqrt(ddot(nvar, grad, 1, grad, 1));
    std::fprintf(stdout, " GRADIENT NORM   %12.6f\n", gnorm);
    cosine = cosine * one;
    std::fprintf(stdout, " DIRECTION COSINE%12.6f\n", cosine);
    if ((iloop - 1) % 50 == 0)
      to_screen(" Structure  No.  Distance from A - B  Heat of Formation Gradient Norm  Cosine of angle");
    if (iflag == 1) {
      char buf[200];
      std::snprintf(buf, sizeof(buf), "  SECOND %5d%18.6f%18.6f%16.6f%16.6f%16.6f",
                    iloop, sum, funct1, gnorm, cosine, 0.0);
      line = buf;
    } else {
      char buf[200];
      std::snprintf(buf, sizeof(buf), "  FIRST  %5d%18.6f%18.6f%16.6f%16.6f%16.6f",
                    iloop, sum, funct1, gnorm, cosine, 0.0);
      line = buf;
    }
    to_screen(trim(line));
    xyz.assign(3 * numat, 0.0);
    for (i = 1; i <= numat; ++i)
      for (k = 1; k <= 3; ++k) xyz[(i - 1) * 3 + (k - 1)] = geo[k][i];
    {
      std::vector<double> geo_buf(3 * numat, 0.0);
      xyzint(xyz.data(), numat, na.data(), nb.data(), nc.data(), 1.0, geo_buf.data());
      for (i = 1; i <= numat; ++i)
        for (k = 1; k <= 3; ++k) geo[k][i] = geo_buf[(i - 1) * 3 + (k - 1)];
    }
    geout(6);
    for (i = 0; i < (int)na.size(); ++i) na[i] = 0;
    for (i = 1; i <= numat; ++i)
      for (k = 1; k <= 3; ++k) geo[k][i] = xyz[(i - 1) * 3 + (k - 1)];
    // F90: if (cosine < 0.0) i = 0  (dead assignment, kept verbatim)
    if (cosine < 0.0) i = 0;
    if (swap > 2.9 || (iloop > 3 && cosine < 0.0) || escf > eold) {
      if (swap > 2.9) {
        swap = 0.0;
      } else {
        swap = 0.5;
      }
      // SWAP REACTANT AND PRODUCT AROUND
      bool finish = (gok[0] && gok[1] && cosine < -1.8);
      if (finish) {
        mopend("BOTH SYSTEMS ARE ON THE SAME SIDE OF THE TRANSITION STATE.");
        return;
      }
      time2 = second(2);
      std::fprintf(stdout, " TIME=%9.2f\n", time2 - time1);
      time1 = time2;
      std::fprintf(stdout, "  REACTANTS AND PRODUCTS SWAPPED AROUND\n");
      iflag = 1 - iflag;
      one = -1.0;
      eold = escf;
      if (gnorm > 10.0) gok[iflag] = true;
      gnorm = sum;
      for (i = 1; i <= natoms; ++i)
        for (j = 1; j <= 3; ++j) {
          x = geo[j][i];
          geo[j][i] = geoa[j][i];
          geoa[j][i] = x;
        }
      for (i = 1; i <= nvar; ++i) {
        x = xold[i];
        xold[i] = xparam[i];
        xparam[i] = x;
      }
      // SWAP AROUND THE DENSITY MATRICES.
      for (i = 1; i <= linear; ++i) {
        x = pastor[i];
        pastor[i] = pa[i];
        pa[i] = x;
        x = pbstor[i];
        pbstor[i] = pb[i];
        pb[i] = x;
        p[i] = pa[i] + pb[i];
      }
      if (finish) break;
    } else {
      one = 1.0;
    }
  }
  std::fprintf(stdout, " AT END OF REACTION\n");
  gold = std::sqrt(ddot(nvar, grad, 1, grad, 1));
  compfg(xparam, true, funct1, true, grad, true);
  if (moperr) return;
  gnorm = std::sqrt(ddot(nvar, grad, 1, grad, 1));
  for (i = 1; i <= nvar; ++i) grold[i] = xparam[i];
  writmo();
  if (moperr) return;

  // THE GEOMETRIES HAVE (A) BEEN OPTIMIZED CORRECTLY, OR
  //                     (B) BOTH ENDED UP ON THE SAME SIDE OF THE T.S.
  // TRANSITION STATE LIES BETWEEN THE TWO GEOMETRIES
  c1 = gold / (gold + gnorm);
  c2 = 1.0 - c1;
  for (i = 1; i <= nvar; ++i)
    xparam[i] = c1 * grold[i] + c2 * xold[i];
  step = 0.0;
  sum = gnorm;
  for (i = 1; i <= nvar; ++i) grold[i] = grad[i];
  compfg(xparam, true, funct1, true, grad, true);
  if (moperr) return;
  gnorm = std::sqrt(ddot(nvar, grad, 1, grad, 1));
  if (gnorm < sum) {
    // Gradient has been minimized - write out "better" ts
    std::fprintf(stdout, " BEST ESTIMATE GEOMETRY OF THE TRANSITION STATE\n");
    std::fprintf(stdout, "\n\n          C1=%8.3f     C2=%8.3f\n", c1, c2);
    escf = funct1;
    writmo();
  } else {
    // Gradient is worse - restore last point
    std::fill(grad.begin(), grad.end(), gold);
    for (i = 1; i <= nvar; ++i) xparam[i] = grold[i];
  }
  iflepo = -1;  // Prevent any more calls to writmo
}

void dock(std::vector<std::vector<double>>& geoa,
          std::vector<std::vector<double>>& geo, double& dist) {
  using namespace molkst_C;

  int i, j, jr, k, l;
  double ca, sa, sum, summ, x;
  static double sum_a[3] = {0.0, 0.0, 0.0};
  static int irot[3][4] = {{0, 0, 0, 0}, {0, 1, 1, 2}, {0, 2, 3, 3}};  // 1-based [1..2][1..3]
  static bool first = true;

  if (first) {
    first = false;
    for (j = 1; j <= numat - id; ++j)
      for (k = 1; k <= 3; ++k) sum_a[k - 1] = sum_a[k - 1] + geo[k][j];
    for (k = 1; k <= 3; ++k) sum_a[k - 1] = sum_a[k - 1] / numat;
  }
  for (j = 1; j <= numat - id; ++j)
    for (k = 1; k <= 3; ++k) geo[k][j] = geo[k][j] - sum_a[k - 1];
  double sum_b[3] = {0.0, 0.0, 0.0};
  for (j = 1; j <= numat - id; ++j)
    for (k = 1; k <= 3; ++k) sum_b[k - 1] = sum_b[k - 1] + geoa[k][j];
  sum = 0.0;
  for (k = 1; k <= 3; ++k) sum_b[k - 1] = sum_b[k - 1] / numat;
  for (j = 1; j <= numat - id; ++j) {
    for (k = 1; k <= 3; ++k) geoa[k][j] = geoa[k][j] - sum_b[k - 1];
    sum = sum + (geo[1][j] - geoa[1][j]) * (geo[1][j] - geoa[1][j]) +
                (geo[2][j] - geoa[2][j]) * (geo[2][j] - geoa[2][j]) +
                (geo[3][j] - geoa[3][j]) * (geo[3][j] - geoa[3][j]);
  }
  for (l = 3; l >= -3; --l) {
    // DOCKING IS DONE IN STEPS OF 16, 4, AND 1 DEGREES AT A TIME.
    ca = std::cos(std::pow(4.0, l - 1) * 0.0174532925199432957);
    sa = std::sqrt(std::fabs(1.0 - ca * ca));
    for (j = 1; j <= 3; ++j) {
      int ir = irot[1][j];
      jr = irot[2][j];
      for (i = 1; i <= 40; ++i) {
        summ = 0.0;
        for (k = 1; k <= numat; ++k) {
          x = ca * geoa[ir][k] + sa * geoa[jr][k];
          geoa[jr][k] = -sa * geoa[ir][k] + ca * geoa[jr][k];
          geoa[ir][k] = x;
          summ = summ + (geo[1][k] - geoa[1][k]) * (geo[1][k] - geoa[1][k]) +
                        (geo[2][k] - geoa[2][k]) * (geo[2][k] - geoa[2][k]) +
                        (geo[3][k] - geoa[3][k]) * (geo[3][k] - geoa[3][k]);
        }
        if (summ > sum) {
          if (i > 1) {
            // Fortran: go to 1000 — undo last rotation, try opposite direction
            sa = -sa;
            for (k = 1; k <= numat; ++k) {
              x = ca * geoa[ir][k] + sa * geoa[jr][k];
              geoa[jr][k] = -sa * geoa[ir][k] + ca * geoa[jr][k];
              geoa[ir][k] = x;
            }
            break;
          }
          sa = -sa;
        } else {
          sum = summ;
        }
        if (i == 1) sum = summ;
      }
    }
  }
  for (j = 1; j <= numat - id; ++j) {
    for (k = 1; k <= 3; ++k) geo[k][j] = geo[k][j] + sum_a[k - 1];
    for (k = 1; k <= 3; ++k) geoa[k][j] = geoa[k][j] + sum_a[k - 1];
  }
  dist = std::sqrt(sum);
}
