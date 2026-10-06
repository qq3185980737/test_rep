// dipind.cpp — C++ translation of MOPAC 2016 "dipind" (dipole.F90-derived
// subroutine originally resident in static_polarizability.F90, split into its
// own translation unit to keep the dependency closure small).
#include "dipind.h"

#include <cmath>
#include <vector>

#include "chrge.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "gmetry.h"
#include "ijbo.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

void dipind(std::vector<double>& dipvec) {
  // MODIFICATION OF DIPOLE SUBROUTINE FOR USE IN THE CALCULATION OF THE
  // INDUCED DIPOLES FOR POLARIZABILITIES.
  std::vector<double> q2(numat + 1, 0.0);
  chrge(p, q2);
  std::vector<double> q(numat + 1, 0.0);
  for (int i = 1; i <= numat; ++i) q[i] = tore[nat[i]] - q2[i];

  std::vector<std::vector<double>> coord(4, std::vector<double>(numat + 1, 0.0));
  gmetry(geo, coord);

  static bool first = true;
  static bool chargd = false;
  static double wtmol = 0.0;
  if (first) {
    wtmol = 0.0;
    double sum = 0.0;
    for (int i = 1; i <= numat; ++i) {
      wtmol = wtmol + ams[nat[i]];
      sum = sum + q[i];
    }
    chargd = std::fabs(sum) > 0.5;
    first = false;
  }
  if (chargd) {
    // NEED TO RESET ION'S POSITION SO THAT THE CENTER OF MASS IS AT THE ORIGIN.
    std::vector<double> center(4, 0.0);
    for (int i = 1; i <= 3; ++i)
      for (int j = 1; j <= numat; ++j) center[i] = center[i] + ams[nat[j]] * coord[i-1][j];
    for (int i = 1; i <= 3; ++i) center[i] = center[i] / wtmol;
    for (int i = 1; i <= 3; ++i)
      for (int j = 1; j <= numat; ++j) coord[i-1][j] = coord[i-1][j] - center[i];
  }
  std::vector<std::vector<double>> dip(5, std::vector<double>(4, 0.0));
  for (int i = 1; i <= numat; ++i) {
    int ni = nat[i];
    int ia = nfirst[i];
    int l = std::min(nlast[i] - ia, 3);
    double hyfsp = 2.0 * dd[ni] * a0 * fpc_8 * fpc_1 * 1.0e-10;
    if (mozyme) {
      for (int j = 1; j <= l; ++j) {
        int k = ijbo(i, i) + 1 + (j * (j + 1)) / 2;
        dip[j][2] = dip[j][2] - hyfsp * p[k];
      }
    } else {
      for (int j = 1; j <= l; ++j) {
        int k = ((ia + j) * (ia + j - 1)) / 2 + ia;
        dip[j][2] = dip[j][2] - hyfsp * p[k];
      }
    }
    for (int d = 1; d <= 3; ++d) dip[d][1] = dip[d][1] + 4.803 * q[i] * coord[d-1][i];
  }
  for (int d = 1; d <= 3; ++d) dip[d][3] = dip[d][2] + dip[d][1];
  for (int j = 1; j <= 3; ++j) {
    dip[4][j] = std::sqrt(dip[1][j] * dip[1][j] + dip[2][j] * dip[2][j] + dip[3][j] * dip[3][j]);
  }
  dipvec[1] = -dip[1][3];
  dipvec[2] = -dip[2][3];
  dipvec[3] = -dip[3][3];
}
