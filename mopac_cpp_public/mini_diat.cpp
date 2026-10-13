// mini_diat.cpp — diat-only micro test with PM7 parameters (zs=1.260237)
// dhc call 46 exact input: xjuc = (1.54005, -5e-5, -5.000000000021654e-05)
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <vector>
#include "parameters_C.h"
#include "overlaps_C.h"
#include "diat.h"
// local forwarder: calls the real C++ diat driver (no mozy_integral_stubs deps)
static void diat_local(int ni, int nj, double* xjuc, double* smat) {
    std::vector<double> xj(xjuc, xjuc + 3);
    std::vector<std::vector<double>> di(10, std::vector<double>(10, 0.0));
    diat(ni, nj, xj, di);
    for (int i = 1; i <= 9; ++i)
        for (int j = 1; j <= 9; ++j) smat[(j - 1) * 9 + (i - 1)] = di[i][j];
}
int main() {
    using namespace parameters_C;
    natorb[1] = 1; zs[1] = 1.260237; zp[1] = 0.0; zd[1] = 0.0;
    npq[1][1] = 1;
    overlaps_C::cutof1 = 1e4;
    double xjuc[3] = {1.54005, -5e-5, -5.000000000021654e-05};
    double sm[81] = {};
    diat_local(1, 1, xjuc, sm);
    unsigned long long b; std::memcpy(&b, &sm[0], 8);
    printf("di(1,1) = %.17g\nhex     = %016llx\n", sm[0], b);
    return 0;
}
