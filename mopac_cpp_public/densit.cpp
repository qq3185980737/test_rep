// densit.cpp — C++ translation of MOPAC 2016 "densit.F90".
#include "densit.h"
#include <vector>

void densit(std::vector<double>& c1d, int mdim, int norbs, int nocc, double occ,
            int nfract, double fract, std::vector<double>& p, int mode) {
    // F90 densit takes c(norbs,norbs) column-major; C++ core uses row-major [i][j].
    std::vector<std::vector<double>> c(norbs + 1, std::vector<double>(norbs + 1, 0.0));
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= norbs; ++j)
            c[i][j] = c1d[(j - 1) * norbs + (i - 1)];
    densit(c, mdim, norbs, nocc, occ, nfract, fract, p, mode);
}


void densit(const std::vector<std::vector<double>>& c, int, int norbs,
            int nocc, double occ, int nfract, double fract,
            std::vector<double>& p, int mode) {
    int norbs2 = norbs / 2;
    double sign, frac, cnst;
    int nl2, nu2, nl1, nu1;
    if (nocc != 0 && nfract > norbs2 && mode != 2) {
        sign = -1.0; frac = occ - fract; cnst = occ;
        nl2 = nfract + 1; nu2 = norbs;
        nl1 = nocc + 1; nu1 = nfract;
    } else {
        sign = 1.0; frac = fract; cnst = 0.0;
        nl2 = 1; nu2 = nocc;
        nl1 = nocc + 1; nu1 = nfract;
    }
    int l = 0;
    for (int i = 1; i <= norbs; ++i) {
        for (int j = 1; j <= i; ++j) {
            l++;
            double sum2 = 0.0, sum1 = 0.0;
            for (int m = nl2; m <= nu2; ++m) sum2 += c[i][m] * c[j][m];
            sum2 *= occ;
            for (int m = nl1; m <= nu1; ++m) sum1 += c[i][m] * c[j][m];
            p[l] = (sum2 + sum1 * frac) * sign;
        }
        p[l] = cnst + p[l];
    }
}
