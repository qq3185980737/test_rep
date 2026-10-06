// charvi.h — C++ translation of MOPAC 2016 "charvi.F90".
#pragma once
#include <vector>

// charvi: character of a vibrational mode under operation ioper.
double charvi(const std::vector<std::vector<double>>& vects, int jorb,
              int ioper, const std::vector<std::vector<double>>& r, int nvecs);
