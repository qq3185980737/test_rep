// parsav.h — C++ translation.
#pragma once
// n and m are inout (Fortran intent(inout)) — restored from the restart file
// on read (mode==0).
void parsav(int mode, int& n, int& m, double* q, double* r, double* efslst,
            double* xlast, int* iiium);
