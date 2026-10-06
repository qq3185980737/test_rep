// rotate_C.cpp — C++ translation of MOPAC 2016 "rotate_C.F90".
// Fortran uses equivalence ccore(i,1) <-> css_i and ccore(i,2) <-> css2_i;
// here the scalars and the array are kept in sync by callers (ROTATE).
#include "rotate_C.h"

namespace rotate_C {
// thread_local: rotate() runs concurrently under the dcart OpenMP atom-pair
// loop; these Fortran-equivalence scratch globals must be per-thread.
thread_local double css1 = 0.0, csp1 = 0.0, cpps1 = 0.0, cppp1 = 0.0;
thread_local double css2 = 0.0, csp2 = 0.0, cpps2 = 0.0, cppp2 = 0.0;
thread_local double ccore[4][3] = {};  // [4][2] in Fortran; +1 on column for 1-based
}  // namespace rotate_C
