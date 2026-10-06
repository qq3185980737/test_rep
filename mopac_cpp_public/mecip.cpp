// mecip.cpp — correct the total density matrix for the CI effect.
#include "mecip.h"

#include <vector>

#include "common_arrays_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mxm.h"

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;

void mecip() {
    // deltap(i,i) = -2*occa(i); off-diagonals to the left = 0.
    for (int i = 1; i <= nmos; ++i) {
        deltap[i][i] = -occa[i] * 2.0;
        for (int j = 1; j < i; ++j) deltap[i][j] = 0.0;
    }

    for (int id = 1; id <= lab; ++id) {
        for (int jd = 1; jd <= id; ++jd) {
            if (nalmat[id] != nalmat[jd]) continue;
            int ix = 0, iy = 0;
            for (int j = 1; j <= nmos; ++j) {
                ix += std::abs(microa[j][id] - microa[j][jd]);
                iy += std::abs(microb[j][id] - microb[j][jd]);
            }
            if (ix + iy > 2) continue;

            if (ix == 2) {
                int i = 1;
                for (; i <= nmos; ++i)
                    if (microa[i][id] != microa[i][jd]) break;
                int ij = microb[i][id];
                int j = i + 1;
                for (; j <= nmos; ++j) {
                    if (microa[j][id] != microa[j][jd]) break;
                    ij += microa[j][id] + microb[j][id];
                }
                double sum = 0.0;
                for (int k = 1; k <= nstate; ++k)
                    sum += vectci[id + (k - 1) * lab] *
                           vectci[jd + (k - 1) * lab];
                double sgn = (1 - 2 * (ij % 2)) / (double)nstate;
                deltap[j][i] += sum * sgn;
            } else if (iy == 2) {
                int i = 1;
                for (; i <= nmos; ++i)
                    if (microb[i][id] != microb[i][jd]) break;
                int ij = 0;
                int j = i + 1;
                for (; j <= nmos; ++j) {
                    if (microb[j][id] != microb[j][jd]) break;
                    ij += microa[j][id] + microb[j][id];
                }
                ij += microa[j][id];
                double sum = 0.0;
                for (int k = 1; k <= nstate; ++k)
                    sum += vectci[id + (k - 1) * lab] *
                           vectci[jd + (k - 1) * lab];
                double sgn = (1 - 2 * (ij % 2)) / (double)nstate;
                deltap[j][i] += sum * sgn;
            } else {
                double sum = 0.0;
                for (int k = 1; k <= nstate; ++k)
                    sum += vectci[id + (k - 1) * lab] *
                           vectci[id + (k - 1) * lab];
                for (int i = 1; i <= nmos; ++i)
                    deltap[i][i] +=
                        (microa[i][id] + microb[i][id]) * sum / nstate;
            }
        }
    }

    // Symmetrize the upper triangle.
    for (int i = 1; i <= nmos; ++i)
        for (int j = 1; j < i; ++j) deltap[j][i] = deltap[i][j];

    // Back-transform: delta = C(active) * deltap; then P += delta * C^T.
    std::vector<double> delta(norbs * nmos, 0.0);
    // Pack deltap and active MO columns into contiguous column-major buffers.
    std::vector<double> dpack(nmos * nmos, 0.0);
    for (int i = 1; i <= nmos; ++i)
        for (int j = 1; j <= nmos; ++j)
            dpack[(i - 1) * nmos + (j - 1)] = deltap[i][j];
    std::vector<double> cf(norbs * nmos, 0.0);
    for (int col = 1; col <= nmos; ++col)
        for (int row = 1; row <= norbs; ++row)
            cf[(col - 1) * norbs + (row - 1)] = c[row][nelec + col];

    // mxm(c(1,nelec+1), norbs, deltap, nmos, delta, nmos)
    mxm(cf.data(), norbs, dpack.data(), nmos, delta.data(), nmos);

    int ij = 0;
    for (int i = 1; i <= norbs; ++i) {
        for (int j = 1; j <= i; ++j) {
            ++ij;
            double sum = 0.0;
            for (int k = 1; k <= nmos; ++k)
                sum += delta[(k - 1) * norbs + (i - 1)] * c[j][nelec + k];
            p[ij] += sum;
        }
    }
}
