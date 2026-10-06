// aabbcd.cpp — C++ translation of MOPAC 2016 "aabbcd.F90" (1-based faithful).
#include "aabbcd.h"

#include "meci_C.h"

static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}

double aabbcd(const int* iocca1, const int* ioccb1, const int* iocca2,
              const int* ioccb2, int nmos, const double* xy) {
    int i = 1, j, k, l, m, ij;
    for (; i <= nmos; ++i) {
        if (iocca1[i] == iocca2[i]) continue;
        break;
    }
    j = i + 1;
    for (; j <= nmos; ++j) {
        if (iocca1[j] == iocca2[j]) continue;
        break;
    }
    k = 1;
    for (; k <= nmos; ++k) {
        if (ioccb1[k] == ioccb2[k]) continue;
        break;
    }
    l = k + 1;
    for (; l <= nmos; ++l) {
        if (ioccb1[l] == ioccb2[l]) continue;
        break;
    }
    if (i == k && j == l && iocca1[i] != ioccb1[i]) {
        meci_C::ispqr[meci_C::iiloop][meci_C::is] = meci_C::jloop;
        meci_C::is++;
    }
    if (iocca1[i] < iocca2[i]) {
        m = i; i = j; j = m;
    }
    if (ioccb1[k] < ioccb2[k]) {
        m = k; k = l; l = m;
    }
    double xr = xy[xyidx(i, j, k, l, nmos)];
    ij = 1;
    if ((i > k && j > l) || (i <= k && j <= l)) ij = 0;
    if (i > k) ij += iocca1[k] + ioccb1[i];
    if (j > l) ij += iocca2[l] + ioccb2[j];
    if (i > k) { m = i; i = k; k = m; }
    for (int p = i; p <= k; ++p) ij += ioccb1[p] + iocca1[p];
    if (j > l) { m = j; j = l; l = m; }
    for (int p = j; p <= l; ++p) ij += ioccb2[p] + iocca2[p];
    if (ij % 2 == 1) xr = -xr;
    return xr;
}
