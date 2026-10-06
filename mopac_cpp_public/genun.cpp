// genun.cpp — C++ translation of MOPAC 2016 "genun.F90".

#include "genun.h"

#include <algorithm>
#include <cmath>

void genun(std::vector<std::vector<double>>& u, int& n) {
    const double pi = 3.14159265358979323846;
    int nequat = (int)std::sqrt(n * pi);
    int nvert = nequat / 2;
    int nu = 0;
    for (int i = 1; i <= nvert + 1; ++i) {
        double fi = pi * (i - 1) / nvert;
        double z = std::cos(fi), xy = std::sin(fi);
        int nhor = (int)(nequat * xy);
        nhor = std::max(1, nhor);
        for (int j = 1; j <= nhor; ++j) {
            if (nu >= n) break;
            double fj = 2.0 * pi * (j - 1) / nhor;
            ++nu;
            u[1][nu] = std::cos(fj) * xy;
            u[2][nu] = std::sin(fj) * xy;
            u[3][nu] = z;
        }
        if (nu >= n) break;
    }
    n = nu;
}

bool collid(double rw, const std::vector<double>& cw,
            const std::vector<std::vector<double>>& cnbr,
            const std::vector<double>& rnbr, int nnbr, int ishape) {
    if (nnbr <= 0 || ishape == 3) return false;
    for (int i = 1; i <= nnbr; ++i) {
        double sumrad = rw + rnbr[i];
        double v1 = std::fabs(cw[1] - cnbr[1][i]);
        if (v1 >= sumrad) continue;
        double v2 = std::fabs(cw[2] - cnbr[2][i]);
        if (v2 >= sumrad) continue;
        double v3 = std::fabs(cw[3] - cnbr[3][i]);
        if (v3 >= sumrad) continue;
        double sr2 = sumrad * sumrad;
        double dd2 = v1 * v1 + v2 * v2 + v3 * v3;
        if (dd2 < sr2) return true;
    }
    return false;
}
