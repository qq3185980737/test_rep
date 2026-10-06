// static_polarizability.h — C++ translation of "static_polarizability.F90".
#pragma once
#include <vector>

// static_polarizability: sets up and runs the finite-field calculation of
// molecular electric response properties (dipole, polarizability, first and
// second hyperpolarizability) via ffhpol.
void static_polarizability();

// ffhpol: finite-field dipole/polarizability (36 compfg field steps).
void ffhpol();

// dipind: induced dipole from the density matrix (point-charge + one-center
// hybridization terms).
void dipind(std::vector<double>& dipvec);
