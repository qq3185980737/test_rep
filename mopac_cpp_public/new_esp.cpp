// new_esp.cpp — C++ translation of MOPAC 2016 "new_esp.F90".
//
// Computes the electrostatic potential (ESP) on a grid around the molecule
// (nuclear + electronic contributions via STO-6G / Rys quadrature), then fits
// the ESP to point-charge / multipole expansions with a least-squares solve
// (LINPACK SQRDC/SQRSL).  Output is written in Jmol CUBE format or a simple
// .grd grid, depending on the "CUBE" / "ESPGRID" keywords.
#include "new_esp.h"
#include "esp_support.h"
#include "esp_utilities.h"
#include "setupg.h"
#include "sqrdc_sqrsl.h"
#include "mult.h"
#include "reada.h"
#include "molkst_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "esp_C.h"
#include "overlaps_C.h"
#include "common_arrays_C.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

using namespace molkst_C;
using namespace common_arrays_C;
using namespace esp_C;
using namespace funcon_C;
using namespace overlaps_C;
using namespace parameters_C;

namespace {

// Fortran index(s, sub): 1-based position or 0.
int index1(const std::string& s, const std::string& sub) {
  std::size_t p = s.find(sub);
  return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

std::string trim(const std::string& s) {
  std::size_t e = s.find_last_not_of(' ');
  return e == std::string::npos ? std::string() : s.substr(0, e + 1);
}

long nint(double x) { return std::lround(x); }

}  // namespace

void new_esp() {
  bool dubl, sames, cube, espgrid;
  int i, j, ij, nx, ny, nz, iiz, iiy, iix, mini, minj, maxi, maxj, lit, ljt;
  int nshells, k, l, ig, jg, isize, isizes, jsize, jsizes, jgmax, idx, nroots;
  int n, mm, in, jn, pij, kl, lmax, nxyz, nFit = 0, m, el, nUse, nnx, nny, nnz;
  float rtemp, dx, dy, dz;
  float ep[81];

  double x, y, z, delx, dely, delz, iatomx, iatomy, iatomz;
  double jatomx, jatomy, jatomz, px, py, pz, icoeff, jcoeff;
  double dtemp1, dtemp2, dtemp3, rr, iexp, jexp, gamma, gammainverse, kfac;
  double contractionDensity[100], tosp, yz, xyz, rx, u[4], w[4], uu, ww;
  double ttinverse, t, x0, y0, z0, tt, xint, yint, zint;
  double xin[64], yin[64], zin[64], EUpperRange, inorm, jnorm;
  double xmax = -1.0e8, xmin = 1.0e8, ymax = -1.0e8, ymin = 1.0e8;
  double zmax = -1.0e8, zmin = 1.0e8;

  std::vector<int> jpvt, itypes, istart, istop, icoord;
  std::vector<float> aMatrix, qraux, bVector, aVector, rsd, qy, qty, xb;
  std::vector<float> esp_array;
  std::vector<std::vector<double>> norm;
  int ijx[100], ijy[100], ijz[100];
  int itype, jtype;
  std::vector<float> planexy, Scr2;
  std::vector<std::vector<double>> s, vecs;

  int norbs_ = norbs;
  s.assign(norbs_ + 1, std::vector<double>(norbs_ + 1, 0.0));
  vecs.assign(norbs_ + 1, std::vector<double>(norbs_ + 1, 0.0));
  istart.assign(norbs_ + 1, 0);
  istop.assign(norbs_ + 1, 0);
  get_minus_point_five_overlap(s);
  {
    // mult(c, s, vecs, norbs): flatten c, s to mult's row-major layout.
    std::vector<double> c_flat(norbs_ * norbs_, 0.0), s_flat(norbs_ * norbs_, 0.0);
    std::vector<double> v_flat(norbs_ * norbs_, 0.0);
    for (i = 1; i <= norbs_; ++i) {
      for (j = 1; j <= norbs_; ++j) {
        c_flat[(i - 1) * norbs_ + (j - 1)] = c[i][j];
        s_flat[(i - 1) * norbs_ + (j - 1)] = s[i][j];
      }
    }
    mult(c_flat.data(), s_flat.data(), v_flat.data(), norbs_);
    for (i = 1; i <= norbs_; ++i)
      for (j = 1; j <= norbs_; ++j) vecs[i][j] = v_flat[(i - 1) * norbs_ + (j - 1)];
  }

  // For GPU MOPAC: density_for_GPU replaces densit.
  ij = norbs_ * (norbs_ + 1) / 2;
  {
    std::vector<double> v_flat(norbs_ * norbs_, 0.0);
    for (i = 1; i <= norbs_; ++i)
      for (j = 1; j <= norbs_; ++j) v_flat[(i - 1) * norbs_ + (j - 1)] = vecs[i][j];
    density_for_GPU(v_flat.data(), fract, nclose, nopen, 2.0, ij, norbs_, 2,
                    p, 5);
  }
  // F90: call dscal(ij, 2.d0, p, 1); then halve diagonal elements.
  for (i = 1; i <= ij; ++i) p[i] *= 2.0;
  for (i = 1; i <= norbs_; ++i) {
    j = i * (i + 1) / 2;
    p[j] = p[j] / 2.0;
  }

  tosp = 2.0 / std::sqrt(3.141592653589793);

  // Convert to AU and work out the upper and lower bounds of the box.
  for (i = 1; i <= numat; ++i)
    for (j = 1; j <= 3; ++j) coord[j-1][i] = coord[j-1][i] / a0;
  for (i = 1; i <= numat; ++i) {
    if (coord[0][i] > xmax) xmax = coord[0][i];
    if (coord[1][i] > ymax) ymax = coord[1][i];
    if (coord[2][i] > zmax) zmax = coord[2][i];
    if (coord[0][i] < xmin) xmin = coord[0][i];
    if (coord[1][i] < ymin) ymin = coord[1][i];
    if (coord[2][i] < zmin) zmin = coord[2][i];
  }
  xmax = xmax + 4.0; ymax = ymax + 4.0; zmax = zmax + 4.0;
  xmin = xmin - 4.0; ymin = ymin - 4.0; zmin = zmin - 4.0;

  nnx = (int)std::max(15, std::min(25, (int)nint(xmax - xmin)));
  nny = (int)std::max(15, std::min(25, (int)nint(ymax - ymin)));
  nnz = (int)std::max(15, std::min(25, (int)nint(zmax - zmin)));
  planexy.assign(nnx * nny * nnz + 1, 0.0f);
  Scr2.assign(nnx * nny * nnz + 1, 0.0f);
  delx = (xmax - xmin) / (nnx - 1);
  dely = (ymax - ymin) / (nny - 1);
  delz = (zmax - zmin) / (nnz - 1);

  // Set up arrays for handling the STO6G orbitals.
  setupg();
  itypes.assign(norbs_ + 1, 0);
  icoord.assign(norbs_ + 1, 0);
  norm.assign(norbs_ + 1, std::vector<double>(7, 0.0));
  nshells = 0;
  for (i = 1; i <= numat; ++i) {
    nshells = nshells + 1;
    itypes[nshells] = 0;          // "s"-type
    icoord[nshells] = i;
    istart[nshells] = nfirst[i];
    istop[nshells] = nfirst[i];
    for (ig = 1; ig <= 6; ++ig)
      norm[nshells][ig] = std::pow(2.0 * zzz[nshells][ig] / 3.141592653589793, 0.75);
    if (nlast[i] == nfirst[i]) continue;
    nshells = nshells + 1;
    itypes[nshells] = 1;          // "p"-type
    icoord[nshells] = i;
    istart[nshells] = nfirst[i] + 1;
    istop[nshells] = nfirst[i] + 3;
    for (ig = 1; ig <= 6; ++ig)
      norm[nshells][ig] =
          std::pow(128.0 / (3.141592653589793 * 3.141592653589793 * 3.141592653589793),
                   0.25) *
          std::pow(zzz[nshells][ig], 1.25);
  }

  // Contribution to the ESP arising from the nuclei.
  std::fill(planexy.begin(), planexy.end(), 0.0f);
  pij = 0;
  for (iiz = 1; iiz <= nnz; ++iiz) {
    z = zmin + delz * (iiz - 1);
    for (iiy = 1; iiy <= nny; ++iiy) {
      y = ymin + dely * (iiy - 1);
      for (iix = 1; iix <= nnx; ++iix) {
        x = xmin + delx * (iix - 1);
        pij = pij + 1;
        for (i = 1; i <= numat; ++i) {
          dx = (float)(x - coord[0][i]);
          dy = (float)(y - coord[1][i]);
          dz = (float)(z - coord[2][i]);
          planexy[pij] = planexy[pij] +
              (float)(tore[nat[i]] / std::sqrt(std::max(
                  (double)dx * dx + (double)dy * dy + (double)dz * dz, 1.0e-12)));
        }
      }
    }
  }

  // Contribution to the ESP arising from the electrons.
  for (i = 1; i <= nshells; ++i) {
    itype = itypes[i];
    iatomx = coord[0][icoord[i]];
    iatomy = coord[1][icoord[i]];
    iatomz = coord[2][icoord[i]];
    if (itype == 0) {       // "s"-type
      isize = 1; isizes = 1; mini = 1; maxi = 1; lit = 1;
    } else {                // "p"-type
      isize = 3; isizes = 3; mini = 2; maxi = 4; lit = 2;
    }
    for (j = 1; j <= i; ++j) {
      sames = (i == j);
      jtype = itypes[j];
      jatomx = coord[0][icoord[j]];
      jatomy = coord[1][icoord[j]];
      jatomz = coord[2][icoord[j]];
      if (jtype == 0) {     // "s"-type
        jsize = 1; jsizes = 1; minj = 1; maxj = 1; ljt = 1;
      } else {              // "p"-type
        jsize = 3; jsizes = 3; minj = 2; maxj = 4; ljt = 2;
      }
      dtemp1 = iatomx - jatomx;
      dtemp2 = iatomy - jatomy;
      dtemp3 = iatomz - jatomz;
      rr = dtemp1 * dtemp1 + dtemp2 * dtemp2 + dtemp3 * dtemp3;
      ij = 0;
      if (sames) {
        for (k = mini; k <= maxi; ++k)
          for (l = minj; l <= k; ++l) {
            ij = ij + 1;
            ijx[ij] = ixn[k] + jxn[l] + 1;
            ijy[ij] = iyn[k] + jyn[l] + 1;
            ijz[ij] = izn[k] + jzn[l] + 1;
          }
      } else {
        for (k = mini; k <= maxi; ++k)
          for (l = minj; l <= maxj; ++l) {
            ij = ij + 1;
            ijx[ij] = ixn[k] + jxn[l] + 1;
            ijy[ij] = iyn[k] + jyn[l] + 1;
            ijz[ij] = izn[k] + jzn[l] + 1;
          }
      }
      for (ig = 1; ig <= 6; ++ig) {
        iexp = zzz[i][ig];
        icoeff = ccc[i][ig];
        inorm = norm[i][ig];
        jgmax = sames ? ig : 6;
        for (jg = 1; jg <= jgmax; ++jg) {
          dubl = (sames && ig != jg);
          jexp = zzz[j][jg];
          jnorm = norm[j][jg];
          jcoeff = ccc[j][jg];
          gamma = iexp + jexp;
          gammainverse = 1.0 / gamma;
          dtemp1 = iexp * jexp * rr * gammainverse;
          if (dtemp1 < 46.0) {
            kfac = std::exp(-dtemp1);
            px = (iexp * iatomx + jexp * jatomx) * gammainverse;
            py = (iexp * iatomy + jexp * jatomy) * gammainverse;
            pz = (iexp * iatomz + jexp * jatomz) * gammainverse;
            dtemp1 = icoeff * jcoeff * inorm * jnorm * kfac;
            idx = 0;
            if (sames) {
              if (dubl) {
                for (k = 1; k <= isize; ++k)
                  for (l = 1; l <= k; ++l) {
                    idx = idx + 1;
                    contractionDensity[idx] = 2.0 * dtemp1;
                  }
              } else {
                for (k = 1; k <= isize; ++k)
                  for (l = 1; l <= k; ++l) {
                    idx = idx + 1;
                    contractionDensity[idx] = dtemp1;
                  }
              }
            } else {
              for (k = 1; k <= isize; ++k)
                for (l = 1; l <= jsize; ++l) {
                  idx = idx + 1;
                  contractionDensity[idx] = dtemp1;
                }
            }
            for (k = 1; k <= ij; ++k)
              contractionDensity[k] = contractionDensity[k] * tosp * gammainverse;
            pij = 1;
            for (iiz = 1; iiz <= nnz; ++iiz) {
              z = zmin + delz * (iiz - 1);
              dz = (float)(pz - z);
              for (iiy = 1; iiy <= nny; ++iiy) {
                y = ymin + dely * (iiy - 1);
                dy = (float)(py - y);
                yz = (double)dy * dy + (double)dz * dz;
                for (iix = 1; iix <= nnx; ++iix) {
                  x = xmin + delx * (iix - 1);
                  isize = isizes;
                  jsize = jsizes;
                  dx = (float)(px - x);
                  xyz = std::max((double)dx * dx + yz, 1.0e-2);
                  rx = gamma * xyz;
                  nroots = (itype + jtype) / 2 + 1;
                  // Solve the Rys polynomials (J. Rys, M. Dupuis, H.F. King,
                  // J. Comput. Chem. 4, 154 (1983)).
                  rys(rx, nroots, u, w);
                  mm = 0;
                  for (n = 1; n <= nroots; ++n) {
                    uu = gamma * u[n];
                    ww = -w[n];
                    tt = gamma + uu;
                    ttinverse = 1.0 / tt;
                    t = std::sqrt(tt);
                    t = 1.0 / t;
                    x0 = (px * gamma + uu * x) * ttinverse;
                    y0 = (py * gamma + uu * y) * ttinverse;
                    z0 = (pz * gamma + uu * z) * ttinverse;
                    in = mm;
                    for (k = 1; k <= lit; ++k) {
                      for (l = 1; l <= ljt; ++l) {
                        jn = in + l;
                        // Solve the one-dimensional integrals using
                        // Gauss-Hermite quadrature.
                        vint(xint, yint, zint, k, l, x0, y0, z0, iatomx,
                             iatomy, iatomz, jatomx, jatomy, jatomz, t);
                        xin[jn] = xint;
                        yin[jn] = yint;
                        zin[jn] = zint * ww;
                      }
                      in = in + 4;
                    }
                    mm = mm + 16;
                  }
                  for (k = 1; k <= ij; ++k) {
                    nx = ijx[k];
                    ny = ijy[k];
                    nz = ijz[k];
                    dtemp1 = 0.0;
                    mm = 0;
                    for (l = 1; l <= nroots; ++l) {
                      dtemp1 = dtemp1 + xin[nx + mm] * yin[ny + mm] * zin[nz + mm];
                      mm = mm + 16;
                    }
                    ep[k] = (float)(contractionDensity[k] * dtemp1);
                  }
                  kl = 0;
                  for (k = istart[i]; k <= istop[i]; ++k) {
                    lmax = sames ? k : istop[j];
                    for (l = istart[j]; l <= lmax; ++l) {
                      kl = kl + 1;
                      planexy[pij] = planexy[pij] +
                          ep[kl] * (float)p[(k * (k - 1)) / 2 + l];
                    }
                  }
                  pij = pij + 1;
                }
              }
            }
          }
        }
      }
    }
  }
  nxyz = pij - 1;
  Scr2 = planexy;
  EUpperRange = 0.9;
  for (;;) {
    for (i = 1; i <= nxyz; ++i)
      if (std::fabs(Scr2[i]) < EUpperRange && std::fabs(Scr2[i]) > 1.0e-4)
        nFit = nFit + 1;
    if (nFit > nxyz - 10) break;
    EUpperRange = EUpperRange * 2.0;
  }
  EUpperRange = EUpperRange * 0.5;
  pij = 1;
  el = 7 * numat;
  aMatrix.assign(el * el + 1, 0.0f);
  qraux.assign(el + 1, 0.0f);
  bVector.assign(el + 1, 0.0f);
  aVector.assign(el + 1, 0.0f);
  rsd.assign(el + 1, 0.0f);
  qy.assign(el + 1, 0.0f);
  qty.assign(el + 1, 0.0f);
  xb.assign(el + 1, 0.0f);
  jpvt.assign(el + 1, 0);
  jpvt.assign(el + 1, 0);   // all columns free (F90: jpvt = 0)
  std::fill(aMatrix.begin(), aMatrix.end(), 0.0f);
  std::fill(bVector.begin(), bVector.end(), 0.0f);
  i = 0;
  l = 0;
  for (iiz = 1; iiz <= nnz; ++iiz) {
    z = zmin + delz * (iiz - 1);
    for (iiy = 1; iiy <= nny; ++iiy) {
      y = ymin + dely * (iiy - 1);
      for (iix = 1; iix <= nnx; ++iix) {
        x = xmin + delx * (iix - 1);
        l = l + 1;
        if (std::fabs(Scr2[l]) < EUpperRange && std::fabs(Scr2[l]) > 1.0e-4) {
          // evec writes 0-based; shift into 1-based storage.
          evec(aVector, x, y, z, coord, numat);
          for (m = el; m >= 1; --m) aVector[m] = aVector[m - 1];
          for (m = 1; m <= el; ++m) {
            for (n = 1; n <= el; ++n)
              aMatrix[(m - 1) * el + n] =
                  aMatrix[(m - 1) * el + n] + aVector[m] * aVector[n];
            bVector[m] = bVector[m] + aVector[m] * Scr2[l];
          }
          i = i + 1;
        }
        pij = pij + 1;
      }
    }
  }

  // Use the Householder transformation to compute the QR factorization of
  // the aMatrix of size el.
  i = 1;
  sqrdc(aMatrix, el, el, el, qraux, jpvt, aVector, i);
  dtemp2 = std::fabs(aMatrix[1]) * 1.0e-5;
  for (k = 1; k <= el; ++k) {
    dtemp1 = std::fabs(aMatrix[(k - 1) * el + k]);
    if (dtemp1 < dtemp2) break;
  }
  nUse = k - 1;

  // SQRSL applies the output of SQRDC to compute coordinate
  // transformations, projections, and least squares solutions.
  i = 111;
  sqrsl(aMatrix, el, el, nUse, qraux, bVector, qy, qty, aVector, rsd, xb, i,
        j);
  for (j = 1; j <= el; ++j) {
    jpvt[j] = -jpvt[j];
    if (j > nUse) aVector[j] = 0.0f;
  }


  // Untangle the vectors.
  for (j = 1; j <= el; ++j) {
    if (jpvt[j] <= 0) {
      k = -jpvt[j];
      while (k != j) {
        rtemp = aVector[j];
        aVector[j] = aVector[k];
        aVector[k] = rtemp;
        jpvt[k] = -jpvt[k];
        k = jpvt[k];
      }
    }
  }
  xmax = xmax + 4.0; xmin = xmin - 4.0;
  ymax = ymax + 4.0; ymin = ymin - 4.0;
  zmax = zmax + 4.0; zmin = zmin - 4.0;
  i = index1(keywrd, " ESPGR");
  if (i != 0) {
    nnx = (int)std::max(nint(reada(keywrd, i + 8)), (long)60);
    nny = nnx;
    nnz = nnx;
  } else {
    nnx = 60; nny = 60; nnz = 60;
  }
  esp_array.assign(nny * nnz + 1, 0.0f);
  delx = (xmax - xmin) / (nnx - 1);
  dely = (ymax - ymin) / (nny - 1);
  delz = (zmax - zmin) / (nnz - 1);

  cube = (index1(keywrd, " CUBE") != 0);
  espgrid = (index1(keywrd, " ESPGRID") != 0);
  if (cube) {
    // Write out ESP data in format for Jmol.
    std::ofstream grd(trim(jobnam) + ".grd");
    grd << " 4 Density\n";
    grd << " Electron density from Total SCF Density\n";
    char buf[200];
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", numat, xmin,
                  ymin, zmin);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nnx, delx, 0.0,
                  0.0);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nny, 0.0, dely,
                  0.0);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nnz, 0.0, 0.0,
                  delz);
    grd << buf;
    for (i = 1; i <= numat; ++i) {
      std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f%12.6f\n", nat[i],
                    (double)nat[i], coord[0][i], coord[1][i], coord[2][i]);
      grd << buf;
    }
    // Write out electrostatic values.
    for (iix = 1; iix <= nnx; ++iix) {
      x = xmin + delx * (iix - 1);
      pij = 0;
      for (iiy = 1; iiy <= nny; ++iiy) {
        y = ymin + dely * (iiy - 1);
        for (iiz = 1; iiz <= nnz; ++iiz) {
          z = zmin + delz * (iiz - 1);
          pij = pij + 1;
          evec(bVector, x, y, z, coord, numat);
          for (m = el; m >= 1; --m) bVector[m] = bVector[m - 1];
          rtemp = 0.0f;
          for (k = 1; k <= el; ++k) rtemp = rtemp + aVector[k] * bVector[k];
          esp_array[pij] = rtemp;
        }
      }
      for (int q = 1; q <= pij; ++q) {
        std::snprintf(buf, sizeof(buf), "%13.5e", esp_array[q]);
        grd << buf;
        if (q % 5 == 0) grd << "\n";
      }
      if (pij % 5 != 0) grd << "\n";
    }
  } else if (espgrid) {
    std::ofstream grd(trim(jobnam) + ".grd");
    char buf[200];
    std::snprintf(buf, sizeof(buf), "%12.6f%12.6f%12.6f\n", xmin * a0,
                  ymin * a0, zmin * a0);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nnx,
                  delx * a0, delx * a0, delx * a0);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nny,
                  dely * a0, dely * a0, dely * a0);
    grd << buf;
    std::snprintf(buf, sizeof(buf), "%5d%12.6f%12.6f%12.6f\n", nnz,
                  delz * a0, delz * a0, delz * a0);
    grd << buf;
    // Write out electrostatic values.
    for (iiz = 1; iiz <= nnz; ++iiz) {
      z = zmin + delz * (iiz - 1);
      pij = 0;
      for (iiy = 1; iiy <= nny; ++iiy) {
        y = ymin + dely * (iiy - 1);
        for (iix = 1; iix <= nnx; ++iix) {
          x = xmin + delx * (iix - 1);
          pij = pij + 1;
          evec(bVector, x, y, z, coord, numat);
          for (m = el; m >= 1; --m) bVector[m] = bVector[m - 1];
          rtemp = 0.0f;
          for (k = 1; k <= el; ++k) rtemp = rtemp + aVector[k] * bVector[k];
          esp_array[pij] = rtemp;
        }
      }
      for (int q = 1; q <= pij; ++q) {
        std::snprintf(buf, sizeof(buf), "%13.5e\n", esp_array[q]);
        grd << buf;
      }
    }
  }
}
