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
    {
        fprintf(stderr, "[BABBBC] i=%d j=%d ij=%d nmos=%d occa1..2=%+.3f %+.3f ioccb1=%d,%d ioccb2=%d,%d iocca1=%d,%d\n",
                i, j, ij, nmos, meci_C::occa[1], meci_C::occa[2],
                ioccb1[1], ioccb1[2], ioccb2[1], ioccb2[2], iocca1[1], iocca1[2]); fflush(stderr);
        for (int k = 1; k <= nmos; ++k)
            fprintf(stderr, "[BABBBC] k=%d xy_ijkk=%+.6f xy_ikjk=%+.6f t1=%.6f t2=%.6f\n",
                    k, xy[xyidx(i, j, k, k, nmos)], xy[xyidx(i, k, j, k, nmos)],
                    (xy[xyidx(i, j, k, k, nmos)] - xy[xyidx(i, k, j, k, nmos)]) * (ioccb1[k] - meci_C::occa[k]),
                    xy[xyidx(i, j, k, k, nmos)] * (iocca1[k] - meci_C::occa[k])); fflush(stderr);
    }
    if ((ij / 2) * 2 != ij) sum = -sum;
    return sum;
}
