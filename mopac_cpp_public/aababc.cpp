// aababc.cpp — C++ translation of MOPAC 2016 "aababc.F90" (1-based faithful).
#include "aababc.h"

#include "meci_C.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double aababc(const int* iocca1, const int* ioccb1, const int* iocca2,
              int nmos, const double* xy) {
    int i = 1, ij, j;
    for (; i <= nmos; ++i) {
        if (iocca1[i] == iocca2[i]) continue;
        break;
    }
    ij = ioccb1[i];
    j = i + 1;
    for (; j <= nmos; ++j) {
        if (iocca1[j] != iocca2[j]) break;
        ij += iocca1[j] + ioccb1[j];
    }
    double sum = 0.0;
    for (int k = 1; k <= nmos; ++k)
        sum += (xy[xyidx(i, j, k, k, nmos)] - xy[xyidx(i, k, j, k, nmos)]) *
                   (iocca1[k] - meci_C::occa[k]) +
               xy[xyidx(i, j, k, k, nmos)] * (ioccb1[k] - meci_C::occa[k]);
    if (ij % 2 == 1) sum = -sum;
    return sum;
}
