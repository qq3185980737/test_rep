// probe_babbcd.cpp — step trace
#include <cstdio>
static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}
int main() {
    const int NM4 = 4;
    double xy4[4*4*4*4] = {0};
    for (int i = 1; i <= NM4; ++i)
        for (int j = 1; j <= NM4; ++j)
            for (int k = 1; k <= NM4; ++k)
                for (int l = 1; l <= NM4; ++l)
                    xy4[xyidx(i, j, k, l, NM4)] = i + 10.0*j + 100.0*k + 1000.0*l;
    int a1[5] = {0, 0, 0, 0, 0}, b1[5] = {0, 1, 1, 0, 0};
    int a2[5] = {0, 0, 0, 0, 0}, b2[5] = {0, 0, 0, 1, 1};
    const int* iocca1=a1; const int* ioccb1=b1; const int* iocca2=a2; const int* ioccb2=b2;
    int nmos = NM4; const double* xy = xy4;
    int ij = 0;
    int i = 1;
    for (; i <= nmos; ++i) { if (ioccb1[i] >= ioccb2[i]) continue; break; }
    std::printf("i=%d\n", i);
    int j = i + 1;
    for (; j <= nmos; ++j) { if (ioccb1[j] < ioccb2[j]) break; ij += iocca2[j] + ioccb2[j]; }
    std::printf("j=%d ij=%d\n", j, ij);
    ij += iocca2[j];
    std::printf("ij+=a2[j]=%d -> %d\n", iocca2[j], ij);
    int k = 1;
    for (; k <= nmos; ++k) { if (ioccb1[k] <= ioccb2[k]) continue; break; }
    std::printf("k=%d\n", k);
    int l = k + 1;
    for (; l <= nmos; ++l) { if (ioccb1[l] > ioccb2[l]) break; ij += iocca1[l] + ioccb1[l]; }
    std::printf("l=%d ij=%d\n", l, ij);
    ij += iocca1[l];
    std::printf("ij+=a1[l]=%d -> %d\n", iocca1[l], ij);
    double one = ((ij / 2) * 2 == ij) ? 1.0 : -1.0;
    double v = (xy[xyidx(i, k, j, l, nmos)] - xy[xyidx(i, l, k, j, nmos)]) * one;
    std::printf("one=%.1f xy(i,k,j,l)=%.1f xy(i,l,k,j)=%.1f v=%.1f\n", one,
                xy[xyidx(i,k,j,l,nmos)], xy[xyidx(i,l,k,j,nmos)], v);
    return 0;
}
