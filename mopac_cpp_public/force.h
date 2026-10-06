// force.h — C++ translation of MOPAC 2016 "force.F90".
#pragma once
#include <vector>

// phase_lock: force the phase of each eigenvector so the largest
// coefficient is positive. vecs is packed n*n column-major (1-based).
void phase_lock(std::vector<double>& vecs, int n);
void force();
