// freqcy.h — C++ translation of MOPAC 2016 "freqcy.F90".
#pragma once
#include <vector>

// Mass-weight lower-triangle packed Hessian (1-based). fmatrx modified in place.
// oldf (same size) receives the unweighted copy.
void freqcy_mass_weight(std::vector<double>& fmatrx,
                        std::vector<double>& oldf,
                        const std::vector<double>& wtmass, int nvar);
void freqcy(std::vector<double>& fmatrx, std::vector<double>& freq,
            std::vector<double>& travel, bool eorc,
            std::vector<std::vector<double>>& deldip,
            std::vector<double>& ff, std::vector<double>& oldf, bool ts);
