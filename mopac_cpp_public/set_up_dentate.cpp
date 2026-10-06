// set_up_dentate.cpp — C++ translation of "set_up_dentate.F90" (9680 B).
// Contains set_up_dentate(), nsp2_correction(), nsp2_atom_correction(),
// C_triple_bond_C(), Si_O_H_Correction(), Si_O_H_bond_correction().
#include "set_up_dentate.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "bangle.h"
#include "common_arrays_C.h"
#include "extvdw_for_MOZYME.h"
#include "funcon_C.h"
#include "mod_atomradii.h"
#include "molkst_C.h"
#include "MOZYME_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace mod_atomradii;
using namespace molkst_C;
using namespace MOZYME_C;

// ---------------------------------------------------------------------------
void set_up_dentate() {
  // Work out which atoms are connected to which other atoms.  Atoms are
  // assumed to be connected if they are within a certain distance of each
  // other.  On exit: nbonds(i) = number of atoms attached to atom i;
  // ibonds(j,i) = atom numbers of atoms attached to atom i.
  int i, ik, j, jk, k, kl, il, iu, l;
  double rmin, safety;
  double coord1[3];

  nbonds.assign(numat + 1, 0);
  ibonds.assign(16, std::vector<int>(numat + 1, 0));
  radius.assign(numat + 1, 0.0);          // allocate(radius(numat))
  extvdw_for_MOZYME(radius, atom_radius_covalent);

  // Atoms are assumed attached if they are within 1.1 times the sum of
  // their covalent radii.
  for (i = 1; i <= numat; ++i) {
    for (j = 1; j <= i - 1; ++j) {
      rmin = 1.0e6;
      if (id == 0) {
        rmin = (coord[0][i] - coord[0][j]) * (coord[0][i] - coord[0][j]) +
               (coord[1][i] - coord[1][j]) * (coord[1][i] - coord[1][j]) +
               (coord[2][i] - coord[2][j]) * (coord[2][i] - coord[2][j]);
      } else {
        for (ik = -l11; ik <= l11; ++ik)
          for (jk = -l21; jk <= l21; ++jk)
            for (kl = -l31; kl <= l31; ++kl) {
              for (int d = 1; d <= 3; ++d)
                coord1[d - 1] = coord[d - 1][j] + common_arrays_C::tvec[d][1] * ik +
                                common_arrays_C::tvec[d][2] * jk + common_arrays_C::tvec[d][3] * kl;
              double dd = (coord1[0] - coord[0][i]) * (coord1[0] - coord[0][i]) +
                          (coord1[1] - coord[1][i]) * (coord1[1] - coord[1][i]) +
                          (coord1[2] - coord[2][i]) * (coord1[2] - coord[2][i]);
              rmin = std::min(rmin, dd);
            }
      }

      // Apply safety criteria to specific diatomic pairs.
      il = std::min(nat[i], nat[j]);
      iu = std::max(nat[i], nat[j]);
      safety = 1.1;
      switch (il) {
        case 1:
          if (iu == 6) safety = 1.25;  // C-H
          break;
        case 5:
          if (iu == 7) safety = 1.0;   // B-N bonds are unusually short
          break;
        case 6:
          if (iu < 8) safety = 1.2;    // C-C, C-N, C-O
          break;
        case 16:
          if (iu == 16) safety = 1.2;  // S-S
          break;
        default:
          break;
      }
      double rr = safety * (radius[i] + radius[j]);
      if (rmin < rr * rr) {
        if (nbonds[i] < 15 && nbonds[j] < 15) {
          nbonds[i] = nbonds[i] + 1;
          nbonds[j] = nbonds[j] + 1;
          ibonds[nbonds[i]][i] = j;
          ibonds[nbonds[j]][j] = i;
        }
      }
    }
  }

  // Check for H attached to H: a hydrogen with 2+ bonds that are all
  // hydrogens loses all but its bonds to non-hydrogens.
  for (i = 1; i <= numat; ++i) {
    if (nat[i] != 1) continue;
    if (nbonds[i] < 2) continue;
    k = 0;
    for (j = 1; j <= nbonds[i]; ++j) {
      if (nat[ibonds[j][i]] != 1) {
        ++k;
        ibonds[k][i] = ibonds[j][i];
      }
    }
    nbonds[i] = k;
  }
  if (!nijbo.empty()) {
    // Sanity check: don't allow a bond if the atoms are too far apart for
    // nijbo to be positive.
    for (i = 1; i <= numat; ++i) {
      k = 0;
      for (j = 1; j <= nbonds[i]; ++j) {
        l = ibonds[j][i];
        if (nijbo[l][i] > -1) {
          ++k;
          ibonds[k][i] = l;
        }
        nbonds[i] = k;
      }
    }
  }
}

// ---------------------------------------------------------------------------
double nsp2_correction() {
  // Add a molecular mechanics correction to all nitrogen atoms that have
  // exactly three ligands.
  double correction, sum;
  if (!(method_pm6 || method_pm7)) return 0.0;
  correction = 0.0;
  for (int i = 1; i <= numat; ++i) {
    if (nat[i] == 7 && nbonds[i] == 3) {
      int j = 0;
      if (nat[ibonds[1][i]] == 1) j = 1;
      if (nat[ibonds[2][i]] == 1) ++j;
      if (nat[ibonds[3][i]] == 1) ++j;
      if (j < 2) {
        sum = nsp2_atom_correction(coord, i, ibonds[1][i], ibonds[2][i],
                                   ibonds[3][i]);
        correction = correction + sum;
      }
    }
  }
  return correction;
}

// ---------------------------------------------------------------------------
double nsp2_atom_correction(const std::vector<std::vector<double>>& vectors,
                            int n, int i, int j, int k) {
  // Evaluate the penalty for non-planarity - done by working out the three
  // angles about the central atom "n" subtended by the lines to atoms i, j, k.
  double a, b, c, ab, ac, bc, tot, cosa, cosb, cosc;
  a = std::sqrt((vectors[1][n] - vectors[1][i]) * (vectors[1][n] - vectors[1][i]) +
                (vectors[2][n] - vectors[2][i]) * (vectors[2][n] - vectors[2][i]) +
                (vectors[3][n] - vectors[3][i]) * (vectors[3][n] - vectors[3][i]));
  b = std::sqrt((vectors[1][n] - vectors[1][j]) * (vectors[1][n] - vectors[1][j]) +
                (vectors[2][n] - vectors[2][j]) * (vectors[2][n] - vectors[2][j]) +
                (vectors[3][n] - vectors[3][j]) * (vectors[3][n] - vectors[3][j]));
  c = std::sqrt((vectors[1][n] - vectors[1][k]) * (vectors[1][n] - vectors[1][k]) +
                (vectors[2][n] - vectors[2][k]) * (vectors[2][n] - vectors[2][k]) +
                (vectors[3][n] - vectors[3][k]) * (vectors[3][n] - vectors[3][k]));
  ab = std::sqrt((vectors[1][j] - vectors[1][i]) * (vectors[1][j] - vectors[1][i]) +
                 (vectors[2][j] - vectors[2][i]) * (vectors[2][j] - vectors[2][i]) +
                 (vectors[3][j] - vectors[3][i]) * (vectors[3][j] - vectors[3][i]));
  ac = std::sqrt((vectors[1][k] - vectors[1][i]) * (vectors[1][k] - vectors[1][i]) +
                 (vectors[2][k] - vectors[2][i]) * (vectors[2][k] - vectors[2][i]) +
                 (vectors[3][k] - vectors[3][i]) * (vectors[3][k] - vectors[3][i]));
  bc = std::sqrt((vectors[1][j] - vectors[1][k]) * (vectors[1][j] - vectors[1][k]) +
                 (vectors[2][j] - vectors[2][k]) * (vectors[2][j] - vectors[2][k]) +
                 (vectors[3][j] - vectors[3][k]) * (vectors[3][j] - vectors[3][k]));
  // dacos = acos; guard division by zero.
  double denom;
  denom = 2.0 * b * c;
  cosa = denom > 1e-14 ? std::acos(std::max(-1.0, std::min(1.0,
          (b * b + c * c - bc * bc) / denom))) : 0.0;
  denom = 2.0 * a * c;
  cosb = denom > 1e-14 ? std::acos(std::max(-1.0, std::min(1.0,
          (a * a + c * c - ac * ac) / denom))) : 0.0;
  denom = 2.0 * b * a;
  cosc = denom > 1e-14 ? std::acos(std::max(-1.0, std::min(1.0,
          (b * b + a * a - ab * ab) / denom))) : 0.0;
  // tot = difference between the sum of the three angles and 360 degrees,
  // expressed as radians.  4*asin(1) = 2*pi.
  tot = 4.0 * std::asin(1.0) - (cosa + cosb + cosc);
  return -0.5 * std::exp(-10.0 * tot);
}

// ---------------------------------------------------------------------------
double C_triple_bond_C() {
  // Energy contribution from acetylenic bonds - correction to account for
  // the extra stabilization of yne bonds.
  if (!(method_pm6 || method_pm7)) return 0.0;
  int isum = 0;
  for (int i = 1; i <= numat; ++i) {
    if (nat[i] == 6 && nbonds[i] == 2) {
      // Possible triple bond.
      int j = 0;
      double rab = 10.0;
      if (nat[ibonds[1][i]] == 6) {
        j = ibonds[1][i];
        rab = (coord[0][i] - coord[0][j]) * (coord[0][i] - coord[0][j]) +
              (coord[1][i] - coord[1][j]) * (coord[1][i] - coord[1][j]) +
              (coord[2][i] - coord[2][j]) * (coord[2][i] - coord[2][j]);
      }
      if (nat[ibonds[2][i]] == 6) {
        j = ibonds[2][i];
        rab = std::min(rab, (coord[0][i] - coord[0][j]) * (coord[0][i] - coord[0][j]) +
                            (coord[1][i] - coord[1][j]) * (coord[1][i] - coord[1][j]) +
                            (coord[2][i] - coord[2][j]) * (coord[2][i] - coord[2][j]));
      }
      if (j == 0) continue;
      if (rab > 1.65) continue;  // squared distance > 1.65^2? Fortran: rab>1.65
      // Carbon atom i is attached to two other atoms and to carbon atom j,
      // and the i-j distance indicates an acetylenic bond.
      ++isum;
    }
  }
  return isum * 6.0;  // the value "6" was determined empirically
}

// ---------------------------------------------------------------------------
double Si_O_H_Correction() {
  // If an Si-O-H structure, add in a bending perturbation.
  double sum = 0.0;
  for (int i = 1; i <= numat; ++i) {
    if (nat[i] == 8) {
      int O = i, Si = 0, H = 0;
      for (int j = 1; j <= nbonds[i]; ++j) {
        int k = ibonds[j][i];
        if (nat[k] == 14) Si = k;
        if (nat[k] == 1) H = k;
      }
      if (Si != 0 && H != 0)
        sum = sum + Si_O_H_bond_correction(coord, Si, O, H);
    }
  }
  return sum;
}

// ---------------------------------------------------------------------------
double Si_O_H_bond_correction(const std::vector<std::vector<double>>& coord,
                              int Si, int O, int H) {
  // Assuming that an Si-O-H structure has been identified, apply a
  // correction.  The angle should be 115 (reference 125) degrees.
  // Apply a Gaussian when Si-O > 1.7 A and O-H > 1.0 A; Gaussian drops to
  // 0.05 when atoms are no longer considered connected
  // = 2.02 A for Si-O and 1.21 A for O-H.
  double r_Si_O =
      (coord[0][O] - coord[0][Si]) * (coord[0][O] - coord[0][Si]) +
      (coord[1][O] - coord[1][Si]) * (coord[1][O] - coord[1][Si]) +
      (coord[2][O] - coord[2][Si]) * (coord[2][O] - coord[2][Si]) - 1.7 * 1.7;
  double r_O_H =
      (coord[0][O] - coord[0][H]) * (coord[0][O] - coord[0][H]) +
      (coord[1][O] - coord[1][H]) * (coord[1][O] - coord[1][H]) +
      (coord[2][O] - coord[2][H]) * (coord[2][O] - coord[2][H]) - 1.0;
  double ref_angle = 125.0 * pi / 180.0;
  double angle = 0.0;
  bangle(coord, Si, O, H, angle);
  return 15.0 * (angle - ref_angle) * (angle - ref_angle) *
         std::exp(-33.0 * std::max(0.0, r_Si_O)) *
         std::exp(-68.0 * std::max(0.0, r_O_H));
}
