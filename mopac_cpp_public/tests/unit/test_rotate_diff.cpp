// Independent check: rotate() finite-difference (delta=1e-3) vs dhcore's ww output.
// H2O CI2 geometry: O(0,0,0.1173) H(0,±0.7572,-0.4692), O-x displacement.
#include <cstdio>
#include <cmath>
#include <vector>
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "rotate.h"

using namespace common_arrays_C;
using namespace molkst_C;

int main() {
    // minimal setup: mimic what setup does for H2O
    molkst_C::numat = 3;
    molkst_C::norbs = 6;
    molkst_C::nfirst[1] = 1; molkst_C::nlast[1] = 4; molkst_C::nat[1] = 8;
    molkst_C::nfirst[2] = 5; molkst_C::nlast[2] = 5; molkst_C::nat[2] = 1;
    molkst_C::nfirst[3] = 6; molkst_C::nlast[3] = 6; molkst_C::nat[3] = 1;
    double coord[9] = {0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692};
    // O-H1 pair
    double xi[3] = {0, 0, 0.1173}, xj[3] = {0, 0.7572, -0.4692};
    double step = 1e-3;
    double wplus[2027] = {}, wminus[2027] = {};
    double e1b[45] = {}, e2a[45] = {}, de1b[45] = {}, de2a[45] = {};
    double enucp = 0, enucm = 0;
    int kr = 1;
    double xip[3] = {step, 0, 0.1173};
    rotate(8, 1, xip, xj, wplus, kr, e1b, e2a, enucp);
    kr = 1;
    double xim[3] = {-step, 0, 0.1173};
    rotate(8, 1, xim, xj, wminus, kr, de1b, de2a, enucm);
    fprintf(stdout, "rotate-diff O-x OH1 w[1..10]:");
    for (int i = 0; i < 10; ++i) fprintf(stdout, " %+.6f", wplus[i] - wminus[i]);
    fprintf(stdout, "\n");
    fprintf(stdout, "wplus [1..10]:");
    for (int i = 0; i < 10; ++i) fprintf(stdout, " %+.6f", wplus[i]);
    fprintf(stdout, "\nwminus[1..10]:");
    for (int i = 0; i < 10; ++i) fprintf(stdout, " %+.6f", wminus[i]);
    fprintf(stdout, "\nenuc diff: %+.6f\n", enucp - enucm);
    return 0;
}
