// dipind.h — C++ translation of MOPAC 2016 "dipind" (from static_polarizability.F90).
#pragma once
#include <vector>

// dipind: induced dipole from the density matrix (point-charge + one-center
// hybridization terms).
void dipind(std::vector<double>& dipvec);
