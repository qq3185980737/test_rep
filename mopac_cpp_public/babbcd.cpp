// babbcd.cpp — C++ translation of MOPAC 2016 "babbcd.F90".
#include "babbcd.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double babbcd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy) {
    int ij = 0;
    int i = 1;
    for (; i <= nmos; ++i) {
        if (ioccb1[i] >= ioccb2[i]) continue;
        break;
    }
    int j = i + 1;
    for (; j <= nmos; ++j) {
        if (ioccb1[j] < ioccb2[j]) break;
        ij += iocca2[j] + ioccb2[j];
    }
    ij += iocca2[j];
    int k = 1;
    for (; k <= nmos; ++k) {
        if (ioccb1[k] <= ioccb2[k]) continue;
        break;
    }
    int l = k + 1;
    for (; l <= nmos; ++l) {
        if (ioccb1[l] > ioccb2[l]) break;
        ij += iocca1[l] + ioccb1[l];
    }
    ij += iocca1[l];
    double one = ((ij / 2) * 2 == ij) ? 1.0 : -1.0;
    return (xy[xyidx(i, k, j, l, nmos)] - xy[xyidx(i, l, k, j, nmos)]) * one;
}
