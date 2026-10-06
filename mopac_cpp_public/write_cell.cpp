// write_cell.cpp — C++ translation of "write_cell" (moldat.F90 lines 1036-1118).
// Prints unit-cell lengths/angles, empirical formula, density, and
// heat-of-formation/gradient per formula unit.
#include "write_cell.h"

#include <cmath>
#include <cstdio>
#include <string>

#include "common_arrays_C.h"
#include "ef_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "reada.h"
#include "volume.h"

using namespace common_arrays_C;
using namespace ef_C;
using namespace funcon_C;
using namespace molkst_C;

void write_cell(int iprt) {
  static int old_nstep = -1;
  if (iprt < 0 || gui) return;
  if (iprt == 0) {
    if (old_nstep == nstep) return;
    old_nstep = nstep;
  }
  // Unit-cell lengths and angles
  double ta = std::sqrt(tvec[1][1] * tvec[1][1] + tvec[2][1] * tvec[2][1] + tvec[3][1] * tvec[3][1]);
  double tb = std::sqrt(tvec[1][2] * tvec[1][2] + tvec[2][2] * tvec[2][2] + tvec[3][2] * tvec[3][2]);
  double tc = std::sqrt(tvec[1][3] * tvec[1][3] + tvec[2][3] * tvec[2][3] + tvec[3][3] * tvec[3][3]);
  double tab = std::sqrt((tvec[1][1] - tvec[1][2]) * (tvec[1][1] - tvec[1][2]) +
                         (tvec[2][1] - tvec[2][2]) * (tvec[2][1] - tvec[2][2]) +
                         (tvec[3][1] - tvec[3][2]) * (tvec[3][1] - tvec[3][2]));
  double tac = std::sqrt((tvec[1][1] - tvec[1][3]) * (tvec[1][1] - tvec[1][3]) +
                         (tvec[2][1] - tvec[2][3]) * (tvec[2][1] - tvec[2][3]) +
                         (tvec[3][1] - tvec[3][3]) * (tvec[3][1] - tvec[3][3]));
  double tbc = std::sqrt((tvec[1][3] - tvec[1][2]) * (tvec[1][3] - tvec[1][2]) +
                         (tvec[2][3] - tvec[2][2]) * (tvec[2][3] - tvec[2][2]) +
                         (tvec[3][3] - tvec[3][2]) * (tvec[3][3] - tvec[3][2]));
  double talpha = 57.295779513 * std::acos((tb * tb + tc * tc - tbc * tbc) / (2.0 * tc * tb));
  double tbeta = 57.295779513 * std::acos((ta * ta + tc * tc - tac * tac) / (2.0 * ta * tc));
  double tgamma = 57.295779513 * std::acos((ta * ta + tb * tb - tab * tab) / (2.0 * ta * tb));
  tab = 1.0;
  if (keywrd.find(" BCC") != std::string::npos) tab = 2.0;
  // Empirical formula
  int z = 0;
  int i = static_cast<int>(keywrd.find(" Z="));
  if (i != 0) {
    i = static_cast<int>(std::lround(reada(keywrd, i)));
    z = mers[1] * mers[2] * mers[3] * i;
    if (keywrd.find(" BCC") != std::string::npos) z = z / 2;
  } else {
    int nel[101] = {0};
    for (i = 1; i <= numat; ++i) nel[nat[i]] = nel[nat[i]] + 1;
    int j = 0;
    for (i = 1; i <= 100; ++i) {
      if (nel[i] > 0) {
        ++j;
        nel[j] = nel[i];
      }
    }
    int k = 10000;
    for (i = 1; i <= j; ++i) {
      if (nel[i] < k) k = nel[i];
    }
    int m = 0;
    int l = 0;
    for (i = 1; i <= 20; ++i) {
      m = 0;
      for (l = 1; l <= j; ++l) {
        if (std::fabs((i * nel[l]) / (double)k - (i * 1.0 * nel[l]) / (double)k) > 1.e-5) m = 1;
      }
      if (m == 0) break;
    }
    z = k / i;
  }
  double vecs[9];
  for (int ii = 1; ii <= 3; ++ii)
    for (int jj = 1; jj <= 3; ++jj) vecs[(jj - 1) * 3 + (ii - 1)] = tvec[ii][jj];
  double vol = volume(vecs, 3);
  if (mers[1] > 0 && mers[2] > 0 && mers[3] > 0) {
    char buf[320];
    std::snprintf(buf, sizeof(buf),
                  "%4d a, b, c, alpha, beta, gamma:%7.3f%7.3f%7.3f%7.2f%7.2f%7.2f Vol:%8.2f Density:%6.3f HoF:%10.3f Grad:%7.2f",
                  nstep + 1,
                  tab * ta / mers[1], tab * tb / mers[2], tab * tc / mers[3],
                  talpha, tbeta, tgamma,
                  vol / (mers[1] * mers[2] * mers[3]),
                  mol_weight * 1.e24 / fpc_10 / vol,
                  escf / z, gnorm / std::sqrt(z * 1.0));
    line = buf;
    std::printf("%s\n", line.c_str());
  }
}
