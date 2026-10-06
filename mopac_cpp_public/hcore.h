// hcore.h — C++ translation of MOPAC 2016 "hcore.F90".
#pragma once
void hcore();
// wstore: keep or discard a block of two-electron integrals based on cutof2.
// w points into the packed W array; kr is advanced by ilim*ilim when any
// |element| > cutof2, otherwise rolled back by the same amount.
void wstore(double* w, int& kr, int ni, int ilim);
