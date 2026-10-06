// pol_vol.cpp — C++ translation of "pol_vol.F90".
// Extracted from the full polar response module so that the finite-field
// static_polarizability path does not drag in the TDHF machinery.
#include "pol_vol.h"

#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace parameters_C;

double pol_vol(double average) {
  double polarizability = average * a0 * a0 * a0;
  if (polvol[99] > 1e-4) {
    polarizability = polarizability * polvol[99] + polvol[100];
    for (int i = 1; i <= numat; ++i) polarizability += polvol[labels[i]];
  }
  return polarizability / (a0 * a0 * a0);
}
