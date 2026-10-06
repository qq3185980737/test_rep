// mpcsyb.h — C++ translation.
#pragma once
#include <cstdio>
void mpcsyb(double* chr, int kchrge, double eionis, double& dip);

// F90 mpcpop(icok): Mulliken populations; isyb = open SYBYL output stream.
void mpcpop(int icok, FILE* isyb);
