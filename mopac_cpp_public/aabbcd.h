// aabbcd.h — C++ translation of MOPAC 2016 "aabbcd.F90".
// Two microstates differing by two alpha + two beta MOs. Flat column-major
// 4-D xy (nmos^4), F90 1-based semantics (element 0 padding).
// Side effect: under (i==k && j==l && iocca1[i]!=ioccb1[i]) writes
// meci_C::ispqr[iiloop][is] = jloop and increments meci_C::is.
#pragma once
double aabbcd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy);
