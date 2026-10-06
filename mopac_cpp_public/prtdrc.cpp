// prtdrc.cpp — C++ translation of "prtdrc.F90" (19699 B).
#include "prtdrc.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "chrge.h"
#include "common_arrays_C.h"
#include "dot.h"
#include "drc_C.h"
#include "drcout.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "quadr.h"
#include "reada.h"
#include "xyzint.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace drc_C;
using namespace molkst_C;
using namespace parameters_C;

namespace {

int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// Fortran SIGN(a, b): magnitude of a with sign of b.
double sign_f(double a, double b) {
  return b >= 0 ? std::fabs(a) : -std::fabs(a);
}

}  // namespace

void prtdrc(double deltt, std::vector<double>& xparam,
            std::vector<double>& ref, double ekin, double& gtot, double& etot,
            std::vector<double>& velo0,
            const std::vector<std::array<int, 2>>& mcoprt, int ncoprt,
            bool parmax) {
  // SAVE variables (Fortran SAVE across calls).
  static int icalcn = 0, iloop = 0, ione = 0;
  static double refscf = 0.0, totime = 0.0, told2 = 0.0, told1 = 0.0,
                refx = 0.0, tlast = 0.0, tref = 0.0, stept = 0.0, steph = 0.0,
                stepx = 0.0, gtot0 = 0.0, gtot1 = 0.0, escf0 = 0.0, escf1 = 0.0,
                ekin0 = 0.0, ekin1 = 0.0, etot0 = 0.0, etot1 = 0.0,
                xold0 = 0.0, xold1 = 0.0, xold2 = 0.0, xtot0 = 0.0,
                xtot1 = 0.0, xtot2 = 0.0;
  static bool goturn = false, ldrc = false;
  static std::string cotype[3] = {"BL", "BA", "DI"};

  int i, j, k, l, nfract = 0, ii, n;
  std::vector<double> escf3(4, 0.0), ekin3(4, 0.0), charge, xold3(4, 0.0),
      geo, tsteps(200, 0.0), etot3(4, 0.0), xtot3(4, 0.0);
  double sum, deltat, t1, t2, sum1, dh, cc, bb, aa, c1, fract = 0.0, time = 0.0,
      escf2 = 0.0, ekin2 = 0.0, etot2 = 0.0;
  std::string text1(3, ' '), text2(2, ' ');

  // Fortran: allxyz/allvel/xyz3/vel3/allgeo/geo3 are module-level common arrays
  // (allocated at setup).  In C++ they are file-scope statics; guarantee they
  // have storage even when the icalcn==numcal init block below is skipped
  // (first call: both are 0, exactly as in Fortran molkst_C numcal=0).
  if (xyz3.empty()) {
    int sz = 3 * numat;
    vref.assign(sz + 1, 0.0);
    vref0.assign(sz + 1, 0.0);
    allxyz.assign(sz * 3 + 1, 0.0);
    allvel.assign(sz * 3 + 1, 0.0);
    xyz3.assign(sz * 3 + 1, 0.0);
    vel3.assign(sz * 3 + 1, 0.0);
    allgeo.assign(sz * 3 + 1, 0.0);
    geo3.assign(sz * 3 + 1, 0.0);
    parref.assign(sz + 1, 0.0);
  }
  int nvar = molkst_C::nvar;

  if (icalcn != numcal) {
    vref.clear();
    vref0.clear();
    allxyz.clear();
    allvel.clear();
    xyz3.clear();
    vel3.clear();
    allgeo.clear();
    geo3.clear();
    parref.clear();
    int sz = 3 * numat;
    vref.assign(sz + 1, 0.0);
    vref0.assign(sz + 1, 0.0);
    allxyz.assign(sz * 3 + 1, 0.0);
    allvel.assign(sz * 3 + 1, 0.0);
    xyz3.assign(sz * 3 + 1, 0.0);
    vel3.assign(sz * 3 + 1, 0.0);
    allgeo.assign(sz * 3 + 1, 0.0);
    geo3.assign(sz * 3 + 1, 0.0);
    parref.assign(sz + 1, 0.0);
    icalcn = numcal;
    totime = 0.0;
    etot0 = 0.0;
    etot1 = 0.0;
    etot2 = 0.0;
    escf0 = 0.0;
    escf1 = 0.0;
    escf2 = 0.0;
    ekin0 = 0.0;
    ekin1 = 0.0;
    ekin2 = 0.0;
    gtot = 0.0;
    gtot0 = 0.0;
    gtot1 = 0.0;
    refx = 0.0;
    fract = 0.0;
    told2 = 0.0;
    xold0 = 0.0;
    xold1 = 0.0;
    xold2 = 0.0;
    xtot0 = 0.0;
    xtot1 = 0.0;
    xtot2 = 0.0;
    for (i = 1; i <= nvar; ++i) parref[i] = xparam[i];
    etot = escf + ekin;
    iloop = 1;  // Fortran prtdrc.F90 L111
    tlast = 0.0;
    goturn = false;
    sum = 0.0;
    for (i = 1; i <= nvar; ++i) {
      sum = sum + velo0[i] * velo0[i];
      vref0[i] = velo0[i];
      vref[i] = velo0[i];
    }
    ione = 1;
    ldrc = sum > 1.0;
    iloop = 1;
    told1 = 0.0;
    //
    //  DETERMINE TYPE OF PRINT: TIME, ENERGY OR GEOMETRY PRIORITY
    //  OR PRINT ALL POINTS
    //
    stept = 0.0;
    steph = 0.0;
    stepx = 0.0;
    i = index1(keywrd, " T-PRI");
    if (i != 0) {
      //  Check for "=" sign
      j = index1(keywrd.substr(i + 5), " ") + i + 6;
      for (int ii2 = i + 6; ii2 <= j; ++ii2) {
        if (ii2 <= (int)keywrd.size() && keywrd[ii2 - 1] == '=') j = -1;
      }
      if (j < 0)
        stept = reada(keywrd, index1(keywrd, "T-PRIO") + 5);
      else
        stept = 0.1;
      tref = -1.0e-6;
      std::fprintf(stdout,
                   "\n TIME PRIORITY, INTERVAL =%5.2f FEMTOSECONDS\n\n", stept);
    } else if (index1(keywrd, " H-PRI") != 0) {
      i = index1(keywrd, " H-PRI");
      j = index1(keywrd.substr(i + 5), " ") + i + 6;
      for (int ii2 = i + 6; ii2 <= j; ++ii2) {
        if (ii2 <= (int)keywrd.size() && keywrd[ii2 - 1] == '=') j = -1;
      }
      if (j < 0)
        steph = reada(keywrd, index1(keywrd, "H-PRI") + 5);
      else
        steph = 0.1;
      std::fprintf(stdout,
                   "\n KINETIC ENERGY PRIORITY, STEP =%5.2f KCAL/MOLE\n\n",
                   steph);
    } else if (index1(keywrd, " X-PRI") != 0) {
      i = index1(keywrd, " X-PRI");
      j = index1(keywrd.substr(i + 5), " ") + i + 6;
      for (int ii2 = i + 6; ii2 <= j; ++ii2) {
        if (ii2 <= (int)keywrd.size() && keywrd[ii2 - 1] == '=') j = -1;
      }
      if (j < 0)
        stepx = reada(keywrd, index1(keywrd, "X-PRIO") + 5);
      else
        stepx = 0.05;
      std::fprintf(stdout,
                   "\n GEOMETRY PRIORITY, STEP =%7.4f ANGSTROMS\n\n", stepx);
    }
    if (stepx < 1.0e-6 && steph < 1.0e-6 && stept < 1.0e-6) {
      //  Set default: if a DRC, then time-slice, if an IRC then a movement
      //  slice.
      if (index1(keywrd, " DRC") == 0)
        stepx = 0.0;
      else
        stept = 0.1;
    }
    if (index1(keywrd, " RESTART") != 0 && index1(keywrd, "IRC=") == 0) {
      //  Fortran: unformatted read from unit ires of the whole DRC state.
      //  No C++ stream mapping exists for ires yet; kept as a documented
      //  stub.
    }
  }
  if (escf < -1.0e8) {
    //  Fortran: unformatted write of the whole DRC state to unit ires, then
    //  close(ires, status='KEEP') and return.  Stub: no stream mapped.
    return;
  }
  charge.assign(numat + 1, 0.0);
  chrge(p, charge);
  for (i = 1; i <= numat; ++i) charge[i] = tore[nat[i]] - charge[i];
  deltat = deltt * 1.0e15;
  na = na_store;
  geo.assign(3 * numat + 1, 0.0);
  xyzint(xparam.data(), numat, na.data(), nb.data(), nc.data(),
         57.29577951308232, geo.data());
  if (iloop == 1) {
    etot1 = etot0;
    etot0 = etot;
    escf1 = escf;
    escf0 = escf;
    ekin1 = ekin;
    ekin0 = ekin;
    for (j = 1; j <= 3; ++j)
      for (int q = 1; q <= nvar; ++q) {
        allgeo[(q - 1) * 3 + (j - 1)] = geo[q];
        allxyz[(q - 1) * 3 + (j - 1)] = xparam[q];
        allvel[(q - 1) * 3 + (j - 1)] = velo0[q];
      }
  } else {
    for (int q = 1; q <= nvar; ++q) {
      allgeo[(q - 1) * 3 + 2] = allgeo[(q - 1) * 3 + 1];
      allgeo[(q - 1) * 3 + 1] = allgeo[(q - 1) * 3 + 0];
      allgeo[(q - 1) * 3 + 0] = geo[q];
      allxyz[(q - 1) * 3 + 2] = allxyz[(q - 1) * 3 + 1];
      allxyz[(q - 1) * 3 + 1] = allxyz[(q - 1) * 3 + 0];
      allxyz[(q - 1) * 3 + 0] = xparam[q];
      allvel[(q - 1) * 3 + 2] = allvel[(q - 1) * 3 + 1];
      allvel[(q - 1) * 3 + 1] = allvel[(q - 1) * 3 + 0];
      allvel[(q - 1) * 3 + 0] = velo0[q];
    }
  }
  //
  //  FORM QUADRATIC EXPRESSION FOR POSITION AND VELOCITY W.R.T. TIME.
  //
  t1 = std::max(told2, 0.02);
  t2 = std::max(told1, 0.02) + t1;
  for (i = 1; i <= nvar; ++i) {
    quadr(allgeo[(i - 1) * 3 + 2], allgeo[(i - 1) * 3 + 1],
          allgeo[(i - 1) * 3 + 0], t1, t2, geo3[(i - 1) * 3 + 0],
          geo3[(i - 1) * 3 + 1], geo3[(i - 1) * 3 + 2]);
    quadr(allxyz[(i - 1) * 3 + 2], allxyz[(i - 1) * 3 + 1],
          allxyz[(i - 1) * 3 + 0], t1, t2, xyz3[(i - 1) * 3 + 0],
          xyz3[(i - 1) * 3 + 1], xyz3[(i - 1) * 3 + 2]);
    quadr(allvel[(i - 1) * 3 + 2], allvel[(i - 1) * 3 + 1],
          allvel[(i - 1) * 3 + 0], t1, t2, vel3[(i - 1) * 3 + 0],
          vel3[(i - 1) * 3 + 1], vel3[(i - 1) * 3 + 2]);
  }
  etot2 = etot1;
  etot1 = etot0;
  etot0 = etot;
  quadr(etot2, etot1, etot0, t1, t2, etot3[1], etot3[2], etot3[3]);
  ekin2 = ekin1;
  ekin1 = ekin0;
  ekin0 = ekin;
  quadr(ekin2, ekin1, ekin0, t1, t2, ekin3[1], ekin3[2], ekin3[3]);
  escf2 = escf1;
  escf1 = escf0;
  escf0 = escf;
  quadr(escf2, escf1, escf0, t1, t2, escf3[1], escf3[2], escf3[3]);
  gtot1 = gtot0;
  gtot0 = gtot;
  xtot2 = xtot1;
  xtot1 = xtot0;
  xold2 = xold2 + xold1;
  xold1 = xold0;
  //
  //  CALCULATE CHANGE IN GEOMETRY
  //
  xold0 = 0.0;
  l = 0;
  xtot0 = 0.0;
  for (i = 1; i <= numat; ++i) {
    sum = 0.0;
    sum1 = 0.0;
    for (j = 1; j <= 3; ++j) {
      l = l + 1;
      sum1 = sum1 + (allxyz[(l - 1) * 3 + 0] - ref[l]) *
                        (allxyz[(l - 1) * 3 + 0] - ref[l]);
      sum = sum + (allxyz[(l - 1) * 3 + 1] - allxyz[(l - 1) * 3 + 0]) *
                      (allxyz[(l - 1) * 3 + 1] - allxyz[(l - 1) * 3 + 0]);
    }
    //  xtot0 is the change in geometry from the start of the run
    //  xold0 is the change in geometry from the last step
    xold0 = xold0 + std::sqrt(sum);
    xtot0 = xtot0 + std::sqrt(sum1);
  }
  quadr(xtot2, xtot1, xtot0, t1, t2, xtot3[1], xtot3[2], xtot3[3]);
  quadr(xold2, xold2 + xold1, xold2 + xold1 + xold0, t1, t2, xold3[1],
        xold3[2], xold3[3]);
  //**********************************************************************
  //  GO THROUGH THE CRITERIA FOR DECIDING WHETHER OR NOT TO PRINT THIS
  //  POINT.  IF YES, THEN ALSO CALCULATE THE EXACT POINT AS A FRACTION
  //  BETWEEN THE LAST POINT AND THE CURRENT POINT
  //**********************************************************************
  //  NFRACT IS THE NUMBER OF POINTS TO BE PRINTED IN THE CURRENT DOMAIN
  //**********************************************************************
  if (iloop >= 3) {
    fract = -10.0;
    nfract = 1;
    if (std::fabs(steph) > 1.0e-20) {
      //
      //  CRITERION FOR PRINTING RESULTS IS A CHANGE IN HEAT OF FORMATION =
      //  -CHANGE IN KINETIC ENERGY
      //
      if (refscf == 0.0) {
        i = static_cast<int>(escf2 / steph);
        refscf = i * steph;
      }
      if (iloop == 3) refscf = escf1;
      dh = std::fabs(escf1 - refscf);
      if (dh > steph) {
        steph = sign_f(steph, escf1 - refscf);
        nfract = static_cast<int>(std::fabs(dh / steph));
        cc = escf3[1];
        bb = escf3[2];
        aa = escf3[3];
        //  PROGRAMMERS! - BE VERY CAREFUL IF YOU CHANGE THIS FOLLOWING
        //  SECTION.  THERE IS NUMERICAL INSTABILITY IF ABS(BB/AA) IS VERY
        //  LARGE. NEAR INFLECTION POINTS AA CHANGES SIGN.       JJPS
        if (std::fabs(bb / aa) > 30) {
          //  USE LINEAR INTERPOLATION
          for (i = 1; i <= nfract; ++i)
            tsteps[i] = -(cc - (refscf + i * steph)) / bb;
        } else {
          //  USE QUADRATIC INTERPOLATION
          for (i = 1; i <= nfract; ++i) {
            c1 = cc - (refscf + i * steph);
            tsteps[i] = ((-bb) + sign_f(std::sqrt(bb * bb - 4.0 * (aa * c1)),
                                        bb)) /
                        (2.0 * aa);
          }
        }
        fract = -0.1;
        refscf = refscf + nfract * steph;
      }
    } else if (stept != 0.0) {
      //
      //  CRITERION FOR PRINTING RESULTS IS A CHANGE IN TIME.
      //
      if (std::fabs(totime + told2 - tref) > stept) {
        i = static_cast<int>(totime / stept);
        fract = i * stept - totime;
        i = static_cast<int>((told2 + totime) / stept);
        j = static_cast<int>(totime / stept);
        nfract = i - j + ione;
        ione = 0;
        for (i = 1; i <= nfract; ++i) tsteps[i] = fract + i * stept;
        tref = tref + nfract * stept;
      }
    } else if (stepx != 0.0) {
      //
      //  CRITERION FOR PRINTING RESULTS IS A CHANGE IN GEOMETRY.
      //
      //  refx = integral of change in geometry from the start, quantized by
      //  stepx
      if (xold2 + xold1 - refx > stepx) {
        nfract = std::min(200, static_cast<int>((xold2 + xold1 - refx) / stepx));
        cc = xold3[1];
        bb = xold3[2];
        aa = xold3[3];
        c1 = 0.0;  // Fortran leaves c1 undefined here; keep deterministic 0
        sum = bb * bb - 4.0 * (aa * c1);
        if (std::fabs(bb / aa) > 30 || sum < 1.0e-20) {
          //  USE LINEAR INTERPOLATION
          for (i = 1; i <= nfract; ++i)
            tsteps[i] = -(cc - (refx + i * stepx)) / bb;
        } else {
          //  USE QUADRATIC INTERPOLATION
          for (i = 1; i <= nfract; ++i) {
            c1 = cc - (refx + i * stepx);
            tsteps[i] = ((-bb) + sign_f(std::sqrt(sum), bb)) / (2.0 * aa);
          }
        }
        refx = refx + nfract * stepx;
        fract = -0.1;
      }
    } else {
      //
      //  PRINT EVERY POINT.
      //
      fract = 0.0;
    }
    // Fortran prtdrc.F90 L428: "if (fract >= -9.d0)" stays INSIDE the
    // "if (iloop >= 3)" block; the "else if (iloop == 1)" (L516) pairs with
    // the OUTER "if (iloop >= 3)" and prints the initial t=0 row.
    if (fract >= -9.0) {
      std::vector<std::vector<double>> xyz3v(4, std::vector<double>(nvar + 1)),
          geo3v(4, std::vector<double>(nvar + 1)),
          vel3v(4, std::vector<double>(nvar + 1));
      for (int q = 1; q <= nvar; ++q)
        for (int d = 1; d <= 3; ++d) {
          xyz3v[d][q] = xyz3[(q - 1) * 3 + (d - 1)];
          geo3v[d][q] = geo3[(q - 1) * 3 + (d - 1)];
          vel3v[d][q] = vel3[(q - 1) * 3 + (d - 1)];
        }
      if (fract == 0.0 && nfract == 1) {
        text1 = " ";
        text2 = " ";
        ii = 0;
        drcout(xyz3v, geo3v, vel3v, nvar, totime, escf3, ekin3, etot3, xtot3,
               iloop, charge, fract, text1, text2, ii, itemp_1);
        n = 0;
        for (i = 1; i <= ncoprt; ++i) {
          k = mcoprt[i - 1][0];
          j = mcoprt[i - 1][1];
          l = k * 3 - 3 + j;
          if (std::fabs(geo3[(l - 1) * 3 + 2]) > 1.0e-20)
            fract = -geo3[(l - 1) * 3 + 1] / (geo3[(l - 1) * 3 + 2] * 2.0);
          if (fract <= 0.0 || fract >= told2) continue;
          if (geo3[(l - 1) * 3 + 2] > 0.0) text1 = "MIN";
          if (geo3[(l - 1) * 3 + 2] < 0.0) text1 = "MAX";
          text2 = cotype[j - 1];
          if (n == 0) {
            n = n + 1;
            std::fprintf(stdout, "\n%s\n", std::string(80, '*').c_str());
          }
          time = totime + fract;
          drcout(xyz3v, geo3v, vel3v, nvar, time, escf3, ekin3, etot3, xtot3,
                 iloop, charge, fract, text1, text2, k, itemp_1);
        }
        if (n != 0)
          std::fprintf(stdout, "\n%s\n", std::string(80, '*').c_str());
        if (std::fabs(escf3[3]) > 1.0e-20)
          fract = -escf3[3] / (escf3[3] * 2.0);
        if (!goturn && fract > 0.0 && fract < told2 * 1.04 && parmax) {
          goturn = true;
          time = fract + totime;
          if (escf3[3] > 0.0) {
            text1 = "MIN";
            if (ldrc) {
              sum = dot(velo0, vref, nvar) * dot(velo0, vref, nvar) /
                    (dot(velo0, velo0, nvar) * dot(vref, vref, nvar) + 1.0e-10);
              sum1 = dot(velo0, vref0, nvar) * dot(velo0, vref0, nvar) /
                     (dot(velo0, velo0, nvar) * dot(vref0, vref0, nvar) +
                      1.0e-10);
              if (sum1 > 0.1 && std::fabs(sum1 - 1.0) > 1.0e-6)
                std::fprintf(stdout,
                             "\n COEF. OF V(0)            =%8.5f   LAST V(0)%8.5f   HALF-LIFE =%12.5g FEMTOSECOS\n",
                             sum1, sum, (-0.6931472 * time / std::log(sum1)));
            }
            std::fprintf(stdout,
                         "\n\n HALF-CYCLE TIME =%11.3f FEMTOSECONDS\n",
                         time - tlast);
            tlast = time;
            for (int q = 1; q <= nvar; ++q) vref[q] = velo0[q];
          }
          if (escf3[3] < 0.0) text1 = "MAX";
          text2 = " ";
          drcout(xyz3v, geo3v, vel3v, nvar, time, escf3, ekin3, etot3, xtot3,
                 iloop, charge, fract, text1, text2, 0, itemp_1);
        } else {
          goturn = false;
        }
      } else {
        for (i = 1; i <= nfract; ++i) {
          time = totime + tsteps[i];
          text1 = " ";
          text2 = " ";
          drcout(xyz3v, geo3v, vel3v, nvar, time, escf3, ekin3, etot3, xtot3,
                 iloop, charge, tsteps[i], text1, text2, 0, itemp_1);
        }
      }
      n = 0;
      for (i = 1; i <= ncoprt; ++i) {
        k = mcoprt[i - 1][0];
        j = mcoprt[i - 1][1];
        l = k * 3 - 3 + j;
        if (std::fabs(geo3[(l - 1) * 3 + 2]) > 1.0e-20)
          fract = -geo3[(l - 1) * 3 + 1] / (geo3[(l - 1) * 3 + 2] * 2.0);
        if (fract <= 0.0 || fract >= told2) continue;
        if (geo3[(l - 1) * 3 + 2] > 0.0) text1 = "MIN";
        if (geo3[(l - 1) * 3 + 2] < 0.0) text1 = "MAX";
        text2 = cotype[j - 1];
        if (n == 0) {
          n = n + 1;
          std::fprintf(stdout, "\n%s\n", std::string(80, '*').c_str());
        }
        time = totime + fract;
        drcout(xyz3v, geo3v, vel3v, nvar, time, escf3, ekin3, etot3, xtot3,
               iloop, charge, fract, text1, text2, k, itemp_1);
      }
      if (n != 0)
        std::fprintf(stdout, "\n%s\n", std::string(80, '*').c_str());
    }
  }
  // Fortran prtdrc.F90 L516-521: "else if (iloop == 1)" pairs with the
  // OUTER "if (iloop >= 3)" (L328); iloop==1 prints the initial t=0 row
  // unconditionally.
  else if (iloop == 1) {
    text1 = " ";
    text2 = " ";
    time = 0.0;
    std::vector<std::vector<double>> xyz3v(4, std::vector<double>(nvar + 1)),
        geo3v(4, std::vector<double>(nvar + 1)),
        vel3v(4, std::vector<double>(nvar + 1));
    for (int q = 1; q <= nvar; ++q)
      for (int d = 1; d <= 3; ++d) {
        xyz3v[d][q] = xyz3[(q - 1) * 3 + (d - 1)];
        geo3v[d][q] = geo3[(q - 1) * 3 + (d - 1)];
        vel3v[d][q] = vel3[(q - 1) * 3 + (d - 1)];
      }
    drcout(xyz3v, geo3v, vel3v, nvar, time, escf3, ekin3, etot3, xtot3, iloop,
           charge, fract, text1, text2, 0, itemp_1);
  }
  totime = totime + told2;
  told2 = told1;
  told1 = deltat;
  iloop = iloop + 1;
  //  Fortran: endfile(iw); backspace(iw) — no-op for stdout.
  std::fill(na.begin(), na.end(), 0);
}
