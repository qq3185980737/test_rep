// hcore_for_MOZYME.cpp — C++ translation of MOPAC 2016 "hcore_for_MOZYME.F90"
// (415 lines).  Builds the one-electron matrix h (diagonal IPs + external
// field + one-electron integrals from h1elec) and the packed two-electron W
// array for the MOZYME linear-scaling scheme, using ijbo-blocked atom-pair
// loops.  mode: -1 remove moving-atom terms, 0 fresh build, +1 re-add terms.
#include "hcore_for_MOZYME.h"
#include "add_more_interactions.h"
#include "h1elec.h"
#include "rotate.h"
#include "outer1.h"
#include "outer2.h"
#include "solrot.h"
#include "hcore.h"   // declares wstore(double*, int&, int, int)
#include "ijbo.h"
#include "vecprt_for_MOZYME.h"
#include "reada.h"
#include "mopend.h"
#include "molkst_C.h"
#include "cosmo_C.h"
#include "linear_cosmo.h"
#include "chanel_C.h"
#include "funcon_C.h"
#include "overlaps_C.h"
#include "parameters_C.h"
#include "common_arrays_C.h"
#include "MOZYME_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace cosmo_C;
using namespace chanel_C;
using namespace funcon_C;
using namespace overlaps_C;
using namespace parameters_C;
using namespace common_arrays_C;
using namespace MOZYME_C;

namespace {

int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

std::string trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? std::string() : s.substr(0, e + 1);
}

}  // namespace

void hcore_for_MOZYME() {
  const double eps = 1.0e-10;
  static int imol = 0;
  bool calci = false, calcij = false, calcj = false, fldon = false;
  int i1, i2, ii, im1, io1, ione, ired, j, j1, jj, jo1, jred, k, krmax, kro, ks;
  int mm, nj, kr, i, ni, itemp;
  double const_, enuc = 0.0, fldcon, fnuc = 0.0, half, hterme, xf = 0.0, yf = 0.0, zf = 0.0;
  double xj[3];
  double e1b[45], e2a[45];
  double wjd[2025], wkd[2025];
  double di[81], dibits[81];

  bool debug = (index1(keywrd, " HCORE") != 0);
  add_more_interactions();
  if (moperr) return;
  if (imol != numcal) {
    if (debug || id != 0) {
      std::fprintf(stdout, " Overlap Cutoff Distance:%34.6f CUTOFS\n",
                   std::sqrt(cutofs));
      std::fprintf(stdout, " Cutoff for quadrupolar and higher 2-electron integrals:%12.6f CUTOF2\n",
                   std::sqrt(cutof2));
      std::fprintf(stdout, " Cutoff for dipolar integrals:%26.6f CUTOF1\n",
                   std::sqrt(cutof1));
      if (id != 0)
        std::fprintf(stdout, " Madelung cutoff:%39.6f CUTOFP\n", cutofp);
    }
    imol = numcal;
    xf = 0.0; yf = 0.0; zf = 0.0;
    std::string tmpkey = trim(keywrd);
    i = index1(tmpkey, " FIELD(") + index1(tmpkey, " FIELD=(");
    if (i != 0) {
      // Erase all text from tmpkey except FIELD data (Fortran: tmpkey(:i)=" ",
      // tmpkey(itemp:)=" ", commas blanked in place; reada starts at position i).
      tmpkey.replace(0, i, i, ' ');
      itemp = index1(tmpkey, ")");
      if (itemp != 0) tmpkey.replace(itemp - 1, tmpkey.size() - (itemp - 1),
                                     tmpkey.size() - (itemp - 1), ' ');
      // Read in the effective field in X,Y,Z coordinates.
      xf = reada(tmpkey, i);
      i = index1(tmpkey, ",");
      if (i != 0) {
        tmpkey[i - 1] = ' ';
        yf = reada(tmpkey, i);
        i = index1(tmpkey, ",");
        if (i != 0) {
          tmpkey[i - 1] = ' ';
          zf = reada(tmpkey, i);
        }
      }
      std::fprintf(stdout, "\n          THE ELECTRIC FIELD IS %10.5f %10.5f %10.5f VOLTS/ANGSTROM\n\n",
                   xf, yf, zf);
    }
    // const = Ao/(8h)  (h=Hartree = eV/atomic unit)
    const_ = a0 / ev;
    efield[1] = xf * const_;
    efield[2] = yf * const_;
    efield[3] = zf * const_;
  }
  ione = 1;
  if (id != 0) ione = 0;
  if (mode == -1) {
    for (i = 1; i <= mpack; ++i) h[i] = -h[i];
    enuclr = -enuclr;
  } else if (mode == 0) {
    enuclr = 0.0;
  } else {
    for (i = 1; i <= (int)parth.size() - 1; ++i) h[i] = parth[i];
    enuclr = refnuc;
  }
  krmax = n2elec + 101;
  fldon = false;
  if (std::fabs(efield[1]) > eps || std::fabs(efield[2]) > eps ||
      std::fabs(efield[3]) > eps) {
    fldcon = ev / a0;   // = 51.42
    fldon = true;
  }
  kr = 1;
  mm = 0;
  ired = 1;
  if (mode == 0) {
    for (i = 1; i <= mpack; ++i) h[i] = 0.0;
  }
  for (i = 1; i <= numat; ++i) {
    calci = (jopt[ired] == i);
    if (calci && ired < numred) ired = ired + 1;
    ni = nat[i];
    if (mode == 0) {
      // Fill the diagonals, and off-diagonals on the same atom.
      i2 = ijbo(i, i);
      for (i1 = 1; i1 <= iorbs[i]; ++i1) {
        for (j1 = 1; j1 <= i1; ++j1) {
          i2 = i2 + 1;
          h[i2] = 0.0;
          if (fldon) {
            io1 = i1 - 1;
            jo1 = j1 - 1;
            if ((jo1 == 0) && (io1 == 1)) {
              hterme = -a0 * dd[ni] * efield[1] * fldcon;
              h[i2] = hterme;
            }
            if ((jo1 == 0) && (io1 == 2)) {
              hterme = -a0 * dd[ni] * efield[2] * fldcon;
              h[i2] = hterme;
            }
            if ((jo1 == 0) && (io1 == 3)) {
              hterme = -a0 * dd[ni] * efield[3] * fldcon;
              h[i2] = hterme;
            }
          }
        }
        mm = mm + 1;
        h[i2] = uspd[mm];
        if (fldon) {
          fnuc = -(efield[1] * coord[0][i] + efield[2] * coord[1][i] +
                   efield[3] * coord[2][i]) * fldcon;
          h[i2] = h[i2] + fnuc;
        }
      }
    }
    if (fldon) {
      enuclr = enuclr - fnuc * tore[nat[i]];
    }
    // Fill the atom-other atom one-electron matrix <PSI(LAMBDA)|PSI(SIGMA)>.
    jred = 1;
    im1 = i - ione;
    for (j = 1; j <= im1; ++j) {
      half = 1.0;
      if (i == j) half = 0.5;
      calcj = (jopt[jred] == j);
      if (calcj && jred < numred) jred = jred + 1;
      calcij = (calci || calcj || mode == 0);
      nj = nat[j];
      if (id == 0) {
        // Molecular system.
        if (ijbo(i, j) >= 0) {
          if (calcij) {
            h1elec(ni, nj, &coord[0][i], &coord[0][j], di);
            ii = ijbo(i, j);
            if (i == j) {
              for (i1 = 1; i1 <= iorbs[i]; ++i1)
                for (j1 = 1; j1 <= i1; ++j1) {
                  ii = ii + 1;
                  h[ii] = h[ii] + di[(j1 - 1) * 9 + (i1 - 1)];
                }
            } else {
              for (i1 = 1; i1 <= iorbs[i]; ++i1)
                for (j1 = 1; j1 <= iorbs[j]; ++j1) {
                  ii = ii + 1;
                  h[ii] = h[ii] + di[(j1 - 1) * 9 + (i1 - 1)];
                }
            }
            // Calculate the two-electron integrals W; the electron nuclear
            // terms E1B and E2A; and the nuclear-nuclear term ENUC.
            rotate(ni, nj, &coord[0][i], &coord[0][j], &w[kr], kr, e1b, e2a,
                   enuc);
            enuclr = enuclr + enuc;
          } else if (!direct) {
            kr = kr + (natorb[ni] * (natorb[ni] + 1)) / 2 *
                          (natorb[nj] * (natorb[nj] + 1)) / 2;
          }
        } else if (ijbo(i, j) == -2) {
          if (calcij) {
            outer2(ni, nj, &coord[0][i], &coord[0][j], &w[kr], kr, e1b, e2a,
                   enuc, id, semidr);
            enuclr = enuclr + enuc;
          } else if (!semidr) {
            if (natorb[ni] * natorb[nj] > 0) {
              if (natorb[ni] > 1) {
                if (natorb[nj] > 1)
                  kr = kr + 7;
                else
                  kr = kr + 4;
              } else if (natorb[nj] > 1)
                kr = kr + 4;
              else
                kr = kr + 1;
            }
          }
        } else if (calcij) {
          outer1(ni, nj, &coord[0][i], &coord[0][j], &w[kr], kr, e1b, e2a,
                 enuc, 0, semidr);
          enuclr = enuclr + enuc;
        } else if (!semidr) {
          if (natorb[ni] * natorb[nj] > 0) kr = kr + 1;
        }
        if (calcij) {
          // Add on the electron-nuclear attraction term for atom I.
          ii = ijbo(i, i);
          j1 = (iorbs[i] * (iorbs[i] + 1)) / 2;
          for (i1 = 1; i1 <= j1; ++i1) {
            ii = ii + 1;
            h[ii] = h[ii] + e1b[i1 - 1] * half;
          }
          // Add on the electron-nuclear attraction term for atom J.
          ii = ijbo(j, j);
          j1 = (iorbs[j] * (iorbs[j] + 1)) / 2;
          for (i1 = 1; i1 <= j1; ++i1) {
            ii = ii + 1;
            h[ii] = h[ii] + e2a[i1 - 1] * half;
          }
        }
        if (kr > krmax - 100 && i != numat) {
          std::fprintf(stdout, " %d %d\n", kr, krmax);
          std::fprintf(stdout, " Running out of storage for W in HCORE \n");
          std::fprintf(stdout, " NUMBER OF ATOMS CALCULATED FOR 'W': %d\n", i);
          std::fprintf(stdout, " NUMBER OF ATOMS IN SYSTEM: %d\n", numat);
          mopend("Running out of storage for W in HCORE");
          return;
        }
      } else {
        // Solid-state system.
        if (ijbo(i, j) >= 0) {
          if (calcij) {
            for (int q = 0; q < 81; ++q) di[q] = 0.0;
            for (ii = -l1u; ii <= l1u; ++ii)
              for (jj = -l2u; jj <= l2u; ++jj)
                for (k = -l3u; k <= l3u; ++k) {
                  for (int d = 1; d <= 3; ++d)
                    xj[d - 1] = coord[d-1][j] + common_arrays_C::tvec[d][1] * ii + common_arrays_C::tvec[d][2] * jj +
                                common_arrays_C::tvec[d][3] * k;
                  h1elec(ni, nj, &coord[0][i], xj, dibits);
                  for (int q = 0; q < 81; ++q) di[q] += dibits[q];
                }
            ii = ijbo(i, j);
            if (i == j) {
              for (i1 = 1; i1 <= iorbs[i]; ++i1)
                for (j1 = 1; j1 <= i1; ++j1) {
                  ii = ii + 1;
                  h[ii] = h[ii] + di[(j1 - 1) * 9 + (i1 - 1)];
                }
            } else {
              for (i1 = 1; i1 <= iorbs[i]; ++i1)
                for (j1 = 1; j1 <= iorbs[j]; ++j1) {
                  ii = ii + 1;
                  h[ii] = h[ii] + di[(j1 - 1) * 9 + (i1 - 1)];
                }
            }
            kro = kr;
            solrot(ni, nj, &coord[0][i], &coord[0][j], wjd, wkd, kr, e1b, e2a,
                   enuc);
            for (k = kro; k <= kr - 1; ++k) {
              w[k] = wjd[k - kro];
              wk[k] = wkd[k - kro];
            }
            enuclr = enuclr + enuc;
          } else if (natorb[ni] == 1) {
            if (natorb[nj] == 1)
              kr = kr + 1;
            else
              kr = kr + 10;
          } else if (natorb[nj] == 1)
            kr = kr + 10;
          else
            kr = kr + 100;
        } else if (ijbo(i, j) == -2) {
          if (calcij) {
            ks = kr;
            outer2(ni, nj, &coord[0][i], &coord[0][j], &w[kr], kr, e1b, e2a,
                   enuc, id, semidr);
            for (k = ks; k <= kr - 1; ++k) wk[k] = 0.0;
            enuclr = enuclr + enuc;
          } else if (natorb[ni] * natorb[nj] > 0) {
            if (natorb[ni] > 1) {
              if (natorb[nj] > 1)
                kr = kr + 7;
              else
                kr = kr + 4;
            } else if (natorb[nj] > 1)
              kr = kr + 4;
            else
              kr = kr + 1;
          }
        } else if (calcij) {
          wk[kr] = 0.0;
          outer1(ni, nj, &coord[0][i], &coord[0][j], &w[kr], kr, e1b, e2a,
                 enuc, id, semidr);
          enuclr = enuclr + enuc;
        } else if (natorb[ni] * natorb[nj] > 0) {
          kr = kr + 1;
        }
        if (calcij) {
          ii = ijbo(i, i);
          j1 = (iorbs[i] * (iorbs[i] + 1)) / 2;
          for (i1 = 1; i1 <= j1; ++i1) {
            ii = ii + 1;
            h[ii] = h[ii] + e1b[i1 - 1] * half;
          }
          ii = ijbo(j, j);
          j1 = (iorbs[j] * (iorbs[j] + 1)) / 2;
          for (i1 = 1; i1 <= j1; ++i1) {
            ii = ii + 1;
            h[ii] = h[ii] + e2a[i1 - 1] * half;
          }
        }
        if (kr > krmax - 100 && i != numat) {
          std::fprintf(stdout, " %d %d\n", kr, krmax);
          std::fprintf(stdout, " Running out of storage for W in HCORE ");
          std::fprintf(stdout, " NUMBER OF ATOMS CALCULATED FOR 'W': %d\n", i);
          std::fprintf(stdout, " NUMBER OF ATOMS IN SYSTEM: %d\n", numat);
          mopend("Running out of storage for W in HCORE");
        }
      }
    }
    ii = iorbs[i];
    ii = (ii * (ii + 1)) / 2;
    if (id != 0) {
      for (i1 = kr; i1 <= kr + ii * ii - 1; ++i1) wk[i1] = 0.0;
    }
    if (ii != 0) {
      wstore(&w[kr], kr, ni, ii);
    }
  }
  if (mode == -1) {
    for (i = 1; i <= mpack; ++i) parth[i] = -h[i];
    refnuc = -enuclr;
  }
  if (useps) {
    // In the original code the dielectric correction to the core-core
    // interaction is added to ENUCLR via addnucz(phinet, qscnet, qdenet)
    // (linear_cosmo module, closed source).  Stub: no-op placeholder.
    (void)phinet;
    (void)qscnet;
    (void)qdenet;
  }
  kr = kr - 1;
  if (debug) {
    std::fprintf(stdout, "\n\n          ONE-ELECTRON MATRIX FROM HCORE\n");
    if (mode == -1)
      std::fprintf(stdout, "           AFTER REMOVAL OF TERMS FOR MOVING ATOMS\n");
    if (mode == 1)
      std::fprintf(stdout, "           AFTER ADDITION OF TERMS FOR MOVING ATOMS\n");
    vecprt_for_MOZYME(&h[1], norbs);
    if (kr > 2000) {
      std::fprintf(stdout, " THE TWO-ELECTRON MATRIX IS TOO LARGE TO PRINT\n");
    }
    j = std::min(kr, 2000);
    if (id == 0) {
      std::fprintf(stdout, "\n\n          TWO-ELECTRON MATRIX IN HCORE\n\n");
      for (int q = 1; q <= j; ++q) {
        std::fprintf(stdout, "%10.4f", w[q]);
        if (q % 10 == 0) std::fprintf(stdout, "\n");
      }
      if (j % 10 != 0) std::fprintf(stdout, "\n");
    } else {
      std::fprintf(stdout, "\n\n          TWO-ELECTRON J MATRIX IN HCORE\n\n");
      for (int q = 1; q <= j; ++q) {
        std::fprintf(stdout, "%10.4f", w[q]);
        if (q % 10 == 0) std::fprintf(stdout, "\n");
      }
      if (j % 10 != 0) std::fprintf(stdout, "\n");
      std::fprintf(stdout, "\n\n          TWO-ELECTRON K MATRIX IN HCORE\n\n");
      for (int q = 1; q <= j; ++q) {
        std::fprintf(stdout, "%10.4f", wk[q]);
        if (q % 10 == 0) std::fprintf(stdout, "\n");
      }
      if (j % 10 != 0) std::fprintf(stdout, "\n");
    }
  }
}
