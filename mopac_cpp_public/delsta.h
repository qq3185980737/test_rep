// delsta.h — C++ translation of MOPAC 2016 "delsta.F90".
#pragma once
#include <vector>

// Point-charge electrostatic derivative between atoms ii,jj.
void delsta(const std::vector<int>& nat, const std::vector<int>& iorbs,
            const std::vector<double>& p,
            const std::vector<std::vector<double>>& cdi,
            std::vector<double>& dstat, int ii, int jj);
