#pragma once
#include <vector>
constexpr int N_PARS = 32385*5;
namespace copyc6_data { extern const double pars[N_PARS]; }
// c6ab(iat,jat,iadr,jadr,1..3). maxci(iat) = max compressed index per atom.
void copyc6(int maxc, int max_elem,
    std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>& c6ab,
    std::vector<int>& maxci);
