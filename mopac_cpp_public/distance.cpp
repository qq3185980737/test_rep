// distance.cpp — C++ translation of distance/angle/torsion from
// H_bond_correction_bits.F90 (lines 267-299).
#include "distance.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "bangle.h"
#include "common_arrays_C.h"
#include "dihed.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

double distance(int a, int b) {
  double coord1[3];
  if (id == 0) {
    return std::sqrt((coord[0][a] - coord[0][b]) * (coord[0][a] - coord[0][b]) +
                     (coord[1][a] - coord[1][b]) * (coord[1][a] - coord[1][b]) +
                     (coord[2][a] - coord[2][b]) * (coord[2][a] - coord[2][b]));
  }
  double d = 1.0e6;
  for (int ik = -l1u; ik <= l1u; ++ik)
    for (int jk = -l2u; jk <= l2u; ++jk)
      for (int kl = -l3u; kl <= l3u; ++kl) {
        for (int d_ = 0; d_ <= 2; ++d_)
          coord1[d_] = coord[d_][a] + common_arrays_C::tvec[d_ + 1][1] * ik +
                       common_arrays_C::tvec[d_ + 1][2] * jk +
                       common_arrays_C::tvec[d_ + 1][3] * kl;
        d = std::min(d, (coord1[0] - coord[0][b]) * (coord1[0] - coord[0][b]) +
                        (coord1[1] - coord[1][b]) * (coord1[1] - coord[1][b]) +
                        (coord1[2] - coord[2][b]) * (coord1[2] - coord[2][b]));
      }
  return std::sqrt(d);
}

double angle(int a, int b, int c) {
  double ang = 0.0;
  bangle(coord, a, b, c, ang);
  return ang;
}

double torsion(int i, int j, int k, int l) {
  double ang = 0.0;
  dihed(coord, i, j, k, l, ang);
  return ang;
}
