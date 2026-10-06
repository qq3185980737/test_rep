// aabacd.h — C++ translation of MOPAC 2016 "aabacd.F90".
// Two microstates differing by two alpha MOs; xy flat column-major 4-D,
// F90 1-based semantics.
#pragma once
double aabacd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy);
