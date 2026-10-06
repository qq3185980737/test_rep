// rotate_C.h — C++ translation of MOPAC 2016 "rotate_C.F90".
// Fortran: ccore(4,2) with equivalence ccore(i,1) <-> css_i (i=1..4) and
// ccore(i,2) <-> css2_i. Column dimension kept 1-based -> [3].
#pragma once

namespace rotate_C {
// thread_local: rotate() runs concurrently under the dcart OpenMP atom-pair
// loop; these Fortran-equivalence scratch globals must be per-thread.
extern thread_local double css1, csp1, cpps1, cppp1;  // first  interaction (s, sp, pps, ppp)
extern thread_local double css2, csp2, cpps2, cppp2;  // second interaction
extern thread_local double ccore[4][3];
}  // namespace rotate_C
