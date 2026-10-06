// density_for_MOZYME.h — C++ translation of MOPAC 2016.
#pragma once
#include <vector>

// Build packed density p from compressed occupied LMOs.
void density_for_MOZYME(std::vector<double>& p, int mode, int nclose_loc,
                        const std::vector<double>& partpin);
