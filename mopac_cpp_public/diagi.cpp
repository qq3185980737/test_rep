// diagi.cpp
#include "diagi.h"
double diagi(const int* ialpha, const int* ibeta, const double* eiga,
             const double* xy, int nmos) {
    double x = 0.0;
    for (int i = 1; i <= nmos; ++i) {
        if (ialpha[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j <= nmos; ++j) {
            // xy(i,i,j,j) and xy(i,j,i,j) in column-major (nmos,nmos,nmos,nmos)
            size_t ij_jj = (size_t)(j - 1) * nmos * nmos * nmos
                         + (j - 1) * nmos * nmos
                         + (i - 1) * nmos
                         + (i - 1);
            size_t ij_ij = (size_t)(j - 1) * nmos * nmos * nmos
                         + (i - 1) * nmos * nmos
                         + (j - 1) * nmos
                         + (i - 1);
            x += (xy[ij_jj] - xy[ij_ij]) * ialpha[j] * 0.5
                 + xy[ij_jj] * ibeta[j];
        }
    }
    for (int i = 1; i <= nmos; ++i) {
        if (ibeta[i] == 0) continue;
        x += eiga[i];
        for (int j = 1; j <= i - 1; ++j) {
            size_t ij_jj = (size_t)(j - 1) * nmos * nmos * nmos
                         + (j - 1) * nmos * nmos
                         + (i - 1) * nmos
                         + (i - 1);
            size_t ij_ij = (size_t)(j - 1) * nmos * nmos * nmos
                         + (i - 1) * nmos * nmos
                         + (j - 1) * nmos
                         + (i - 1);
            x += (xy[ij_jj] - xy[ij_ij]) * ibeta[j];
        }
    }
    return x;
}
