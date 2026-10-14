// test_disp.cpp
#include <cstdio>
#include <vector>
#include "H_bond_correction_PM6_DH_Dispersion.h"
int main() { std::vector<double> d(10); double v=PM6_DH_Dispersion(false,d); std::printf("disp links=%g PASS\n",v); return 0; }
