// test_gmetry.cpp
#include <cstdio>
#include <vector>
#include "gmetry.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
using namespace molkst_C;
using namespace common_arrays_C;
int main() {
    natoms = 3; numcal = 1; step = 0.0; id = 0; nvar = 0;
    geo.resize(4, std::vector<double>(natoms+1, 0.0));
    geoa.resize(4, std::vector<double>(natoms+1, 0.0));
    coord.resize(4, std::vector<double>(natoms+1, 0.0));
    na.resize(natoms+1, 0); nb.resize(natoms+1, 0); nc.resize(natoms+1, 0);
    labels.resize(natoms+1, 6);
    loc.resize(3, std::vector<int>(nvar+1, 0));
    xparam.resize(nvar+1, 0.0);
    geo[1][1] = 0; geo[2][1] = 0; geo[3][1] = 0;
    na[2] = 1; geo[1][2] = 1.0;
    gmetry(geo, coord);
    std::printf("coord[1][2]=%g PASS\n", coord[1][2]);
    return 0;
}
