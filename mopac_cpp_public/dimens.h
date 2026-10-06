// dimens.h — C++ translation of MOPAC 2016 "dimens.F90".
#pragma once
#include <vector>

// Compute molecular dimensions (skeleton; mutates coord).
void dimens(std::vector<std::vector<double>>& coord, int iw);
