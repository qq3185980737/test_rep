// bldsym.h — C++ translation of MOPAC 2016 "bldsym.F90".
#pragma once

// bldsym: build the 3x3 point-group operation matrix for ioper into
// symmetry_C::elem(.., .., iplace).
void bldsym(int ioper, int iplace);
