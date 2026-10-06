#include "H_bond_correction_bits.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include <cassert>
#include <cmath>
#include <cstdio>
namespace molkst_C { int numat=0; int norbs=0; int id=0; int l1u=0; int l2u=0; int l3u=0; int l11=0; int l21=0; int l31=0; }
namespace funcon_C { double pi=3.14159265358979323846; }
namespace common_arrays_C { std::vector<std::vector<double>> coord(4, std::vector<double>(3,0.0)); std::vector<int> nat(3,0); std::vector<std::vector<double>> tvec; }
int main() {
    common_arrays_C::coord[0][1]=0; common_arrays_C::coord[1][1]=0; common_arrays_C::coord[2][1]=0;
    common_arrays_C::coord[0][2]=3; common_arrays_C::coord[1][2]=4; common_arrays_C::coord[2][2]=0;
    double d = distance_hb(1,2);
    assert(std::fabs(d-5.0)<1e-9);
    printf("distance_hb PASS\n");
    return 0;
}