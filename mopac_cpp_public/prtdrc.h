// prtdrc.h — C++ translation of "prtdrc.F90".
#pragma once
#include <array>
#include <vector>

// prtdrc: prepares to print the geometry, energy and velocity data for
// points in a DRC or IRC calculation, deciding which points fall on the
// print criteria (time / heat-of-formation / geometry priorities).
void prtdrc(double deltt, std::vector<double>& xparam,
            std::vector<double>& ref, double ekin, double& gtot, double etot,
            std::vector<double>& velo0,
            const std::vector<std::array<int, 2>>& mcoprt, int ncoprt,
            bool parmax);
