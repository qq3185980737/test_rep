// diagg2.h — C++ translation of MOPAC 2016 "diagg2.F90".
#pragma once

// Jacobi annihilation of occupied-virtual LMO couplings with adaptive LMO
// expansion.  Pointer-based signature matches addhb.cpp's forward
// declaration: eigv is 0-based over nvirtual eigenvalues; storei/storej are
// (norbs) scratch; iused/latoms are (numat).
void diagg2(int nocc, int nvir, const double* eigv, int* iused,
            char* latoms, int nij, int idiagg, double* storei,
            double* storej);
