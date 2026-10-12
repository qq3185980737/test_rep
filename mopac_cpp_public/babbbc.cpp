// babbbc.cpp — C++ translation of MOPAC 2016 "babbbc.F90".
#include "babbbc.h"

#include "meci_C.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double babbbc(const int* iocca1, const int* ioccb1, const int* ioccb2,
              int nmos, const double* xy) {
    int i = 1;
    for (; i <= nmos; ++i) {
        if (ioccb1[i] == ioccb2[i]) continue;
        break;
    }
    int ij = 0;
    int j = i + 1;
    for (; j <= nmos; ++j) {
        if (ioccb1[j] != ioccb2[j]) break;
        ij += iocca1[j] + ioccb1[j];
    }
    ij += iocca1[j];
    double sum = 0.0;
    for (int k = 1; k <= nmos; ++k)
        sum += (xy[xyidx(i, j, k, k, nmos)] - xy[xyidx(i, k, j, k, nmos)]) *
                   (ioccb1[k] - meci_C::occa[k]) +
               xy[xyidx(i, j, k, k, nmos)] * (iocca1[k] - meci_C::occa[k]);
    if ((ij / 2) * 2 != ij) sum = -sum;
    return sum;
}
