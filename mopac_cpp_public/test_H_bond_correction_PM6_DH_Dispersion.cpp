// test_H_bond_correction_PM6_DH_Dispersion.cpp
#include <cstdio>
#include <vector>
#include "H_bond_correction_PM6_DH_Dispersion.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace common_arrays_C;
int main() {
    numat = 0; id = 0;
    coord.resize(4, std::vector<double>(1, 0.0));
    nat.resize(1, 0); nbonds.resize(1, 0);
    double e = PM6_DH_Disp(0, 0);
    std::printf("PM6_DH_Disp=%g PASS\n", e);
    std::vector<double> d(1, 0.0);
    double e2 = PM6_DH_Dispersion(false, d);
    std::printf("PM6_DH_Dispersion=%g PASS\n", e2);
    return 0;
}
