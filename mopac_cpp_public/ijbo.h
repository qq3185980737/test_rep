// ijbo.h — C++ translation of ijbo.F90 (MOPAC 2016).
#pragma once

// Given atom numbers ii and jj:
//   IJBO = -1 if the atoms are separated by more than SQRT(CUTOF1), else
//   IJBO = -2 if the atoms are separated by more than SQRT(CUTOF2), else
//   IJBO = N, where "N+1" is the address of the starting location in the
//           H, P, and F arrays.
int ijbo(int ii, int jj);