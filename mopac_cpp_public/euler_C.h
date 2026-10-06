// euler_C.h — C++ translation of MOPAC 2016 "euler_C.F90".
//
// Fortran module euler_C holds: id, l1l..l3l, l1u..l3u, tvec(3,3).
// In this C++ port those globals are split by existing architecture:
//   id, l1u, l2u, l3u  -> molkst_C  (molkst_C.h)
//   tvec(3,3)          -> common_arrays_C (vector form, periodic lattice)
//   l1l, l2l, l3l      -> this namespace (the only members defined here)
#pragma once

namespace euler_C {
extern int l1l, l2l, l3l;
}  // namespace euler_C
