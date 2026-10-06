// ciosci.h — C++ translation of MOPAC 2016 "ciosci.F90".
#pragma once
#include <vector>

// Oscillator-strength / dipole matrix between CI states.
// vects(norbs,norbs), oscil(3,lab), conf(lab*lab).
void ciosci(const std::vector<std::vector<double>>& vects, int lroot,
            std::vector<std::vector<double>>& oscil,
            const std::vector<double>& conf_in);
