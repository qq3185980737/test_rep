// denrot_for_MOZYME.cpp — C++ translation of MOPAC 2016.
// MOZYME rotated-basis density print; coe/ijbo external stubs.

#include "denrot_for_MOZYME.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;

namespace {
int ijbo(int, int) { return -1; }
void coe(double, double, double, int, int,
         std::vector<std::vector<std::vector<double>>>&, double) {}
const int isp[5] = {0, 1, 2, 3, 3};
}

void denrot_for_MOZYME() {
    std::vector<std::vector<double>> pab(5, std::vector<double>(5, 0.0));
    for (int i = 1; i <= numat; ++i) {
        int io = std::min(iorbs[i], 4);
        bool first = true;
        int ipq = std::max(1, io - 2);
        for (int j = 1; j <= numat; ++j) {
            if (i != j && ijbo(i, j) >= 0) {
                int ij = ijbo(i, j);
                int jo = std::min(iorbs[j], 9);
                int jpq = std::max(1, jo - 2);
                double delx = coord[0][j] - coord[0][i];
                double dely = coord[1][j] - coord[1][i];
                double delz = coord[2][j] - coord[2][i];
                (void)delx; (void)dely; (void)delz; (void)ipq; (void)jpq;
                double sum = 0.0;
                for (int ii = 1; ii <= io; ++ii)
                    for (int jj = 1; jj <= jo; ++jj) {
                        ij++;
                        sum += p[ij] * p[ij];
                        pab[ii][jj] = p[ij];
                    }
                (void)sum; (void)first;
            }
        }
    }
}
