// dfpsav.h — C++ translation of MOPAC 2016 "dfpsav.F90".
#pragma once
#include <vector>

// Save/restore DFP geometry-optimization state to disk (skeleton).
void dfpsav(double& totime, std::vector<double>& xparam,
            std::vector<double>& gd, std::vector<double>& xlast,
            double& funct1, std::vector<int>& mdfp, std::vector<double>& xdfp);
