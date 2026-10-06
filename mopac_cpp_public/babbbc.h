// babbbc.h — C++ translation of MOPAC 2016 "babbbc.F90".
// Two microstates differing by one beta electron; xy flat column-major
// 4-D (nmos^4), F90 1-based semantics.
#pragma once
double babbbc(const int* iocca1, const int* ioccb1, const int* ioccb2,
              int nmos, const double* xy);
