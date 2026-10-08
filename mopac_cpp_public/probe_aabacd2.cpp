// probe_aabacd2.cpp — step trace of aabacd logic
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
    const int* iocca1; const int* ioccb1; const int* iocca2; const int* ioccb2;
    int a1[5] = {0, 1, 1, 0, 0}, a2[5] = {0, 0, 0, 1, 1};
    int b1[5] = {0, 0, 0, 0, 0}, b2[5] = {0, 0, 0, 0, 0};
    iocca1=a1; ioccb1=b1; iocca2=a2; ioccb2=b2;
    int nmos = NM4; const double* xy = xy4;
    int ij = 0, i = 1, j, k, l;
    for (; i <= nmos; ++i) { if (iocca1[i] >= iocca2[i]) continue; break; }
    std::printf("i=%d\n", i);
    j = i + 1;
    for (; j <= nmos; ++j) { if (iocca1[j] < iocca2[j]) break; ij += iocca2[j] + ioccb2[j]; }
    std::printf("j=%d ij=%d\n", j, ij);
    k = 1;
    for (; k <= nmos; ++k) { if (iocca1[k] <= iocca2[k]) continue; break; }
    std::printf("k=%d\n", k);
    l = k + 1;
    for (; l <= nmos; ++l) { if (iocca1[l] > iocca2[l]) break; ij += iocca1[l] + ioccb1[l]; }
    std::printf("l=%d ij=%d\n", l, ij);
    ij += ioccb2[i] + ioccb1[k];
    std::printf("ij_final=%d xy(i,k,j,l)=%.1f xy(i,l,k,j)=%.1f\n", ij,
                xy[xyidx(i,k,j,l,nmos)], xy[xyidx(i,l,k,j,nmos)]);
    double sum = xy[xyidx(i, k, j, l, nmos)] - xy[xyidx(i, l, k, j, nmos)];
    if (ij % 2 == 1) sum = -sum;
    std::printf("sum=%.3f\n", sum);
    return 0;
}
