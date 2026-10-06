// deri0.cpp — C++ translation of MOPAC 2016 "deri0.F90".
#include "deri0.h"

#include <algorithm>
#include <cmath>

void deri0(const std::vector<double>& e, int n,
           std::vector<double>& scalar, std::vector<double>& diag,
           double fract, const std::vector<int>& nbo) {
    const double shift = 2.36;
    int nopen = nbo[1] + nbo[2];
    const double cnst = 1e-3;
    int l = 1;
    if (nbo[2] > 0 && nbo[1] > 0) {
        for (int j = 1; j <= nbo[1]; ++j) {
            if (nopen - nbo[1] > 0) {
                for (int m = nbo[1] + 1; m <= nopen; ++m)
                    diag[l + (m - (nbo[1] + 1))] = (e[m] - e[j]) / (2.0 - fract + cnst);
                l = nopen - nbo[1] + l;
            }
        }
    }
    if (nbo[3] > 0 && nbo[1] > 0) {
        for (int j = 1; j <= nbo[1]; ++j) {
            if (n - nopen > 0) {
                for (int m = nopen + 1; m <= n; ++m)
                    diag[l + (m - (nopen + 1))] = (e[m] - e[j]) / 2.0;
                l = n - nopen + l;
            }
        }
    }
    if (nbo[3] != 0 && nbo[2] != 0) {
        for (int j = nbo[1] + 1; j <= nopen; ++j) {
            if (n - nopen > 0) {
                for (int m = nopen + 1; m <= n; ++m)
                    diag[l + (m - (nopen + 1))] = (e[m] - e[j]) / (fract + cnst);
                l = n - nopen + l;
            }
        }
    }
    for (int i = 1; i <= l - 1; ++i)
        scalar[i] = std::sqrt(1.0 / std::max(0.3 * diag[i], diag[i] - shift));
}
