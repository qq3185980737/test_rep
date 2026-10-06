// overlaps_C.cpp — storage for Fortran module "overlaps_C".
// Declarations follow overlaps_C.h (1-based indexing).
#include "overlaps_C.h"

namespace overlaps_C {

double ccc[N_AO + 1][7] = {};
double zzz[N_AO + 1][7] = {};
double allz[7][7][3] = {};
double allc[7][7][3] = {};
double fact[18] = {1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
                   362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0,
                   87178291200.0, 1307674368000.0, 20922789888000.0,
                   355687428096000.0};

// Overlap / two-electron integral cutoffs (defaults used by ijbo & diagg).
double cutof1 = 1.0e-2;
double cutof2 = 1.0e-4;

}  // namespace overlaps_C
