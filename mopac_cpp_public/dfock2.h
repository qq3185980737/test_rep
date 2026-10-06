// dfock2.h — C++ translation of MOPAC 2016 "dfock2.F90".
#pragma once
// Adds 2-electron 2-center repulsion contribution to the Fock matrix
// derivative (NDDO). All packed arrays 1-based (index 0 padding).
void dfock2(double* f, const double* ptot, const double* p, const double* w,
            int numat, const int* nfirst, const int* nlast, int nati);
