// babbcd.h — C++ translation of MOPAC 2016 "babbcd.F90".
// Two microstates differing by two beta MOs; xy is a flat column-major
// 4-D array (nmos^4) with F90 1-based semantics (element 0 padding).
#pragma once
double babbcd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy);
