// ccrep.h — C++ translation of MOPAC 2016 "ccrep.F90".
#pragma once

// ccrep: core-core repulsion between atoms ni,nj.
// Input r (Angstrom) is converted; gab is the monopole factor.
// enuclr is the output repulsion.
void ccrep(int ni, int nj, double& r, double gab, double& enuclr);
