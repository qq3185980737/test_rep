// probe_aabacd.cpp
#include <cstdio>
#include "aabacd.h"
static inline int xyidx(int a, int b, int c, int d, int nmos) {
    return (d - 1) * nmos * nmos * nmos + (c - 1) * nmos * nmos +
           (b - 1) * nmos + (a - 1);
}
int main() {
    const int NM4 = 4;
    double xy4[4*4*4*4];
    for (int i = 0; i < 256; ++i) xy4[i] = 0;
    for (int i = 1; i <= NM4; ++i)
        for (int j = 1; j <= NM4; ++j)
            for (int k = 1; k <= NM4; ++k)
                for (int l = 1; l <= NM4; ++l)
                    xy4[xyidx(i, j, k, l, NM4)] = i + 10.0*j + 100.0*k + 1000.0*l;
    int a1[5] = {0, 1, 1, 0, 0}, a2[5] = {0, 0, 0, 1, 1};
    int b1[5] = {0, 0, 0, 0, 0}, b2[5] = {0, 0, 0, 0, 0};
    double v = aabacd(a1, b1, a2, b2, NM4, xy4);
    double exp = xy4[xyidx(3, 1, 4, 2, NM4)] - xy4[xyidx(3, 2, 4, 1, NM4)];
    std::printf("v=%.3f exp=%.3f xy(3,1,4,2)=%.1f xy(3,2,4,1)=%.1f\n", v, exp,
                xy4[xyidx(3,1,4,2,NM4)], xy4[xyidx(3,2,4,1,NM4)]);
    return 0;
}
