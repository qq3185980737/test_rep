// aababc.h — C++ translation of MOPAC 2016 "aababc.F90".
// Two microstates differing by one alpha electron; xy flat column-major
// 4-D, F90 1-based semantics.
#pragma once
double aababc(const int* iocca1, const int* ioccb1, const int* iocca2,
              int nmos, const double* xy);
