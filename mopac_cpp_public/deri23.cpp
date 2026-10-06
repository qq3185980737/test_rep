// deri23.cpp — C++ translation of MOPAC 2016 "deri23.F90".
#include "deri23.h"

#include <algorithm>
#include <cmath>

#include "meci_C.h"
#include "molkst_C.h"

using namespace meci_C;
using namespace molkst_C;

void deri23(const std::vector<double>& f,
            const std::vector<double>& fd,
            const std::vector<double>& e,
            const std::vector<double>& fci,
            std::vector<std::vector<double>>& cmo,
            std::vector<double>& emo, int minear, int ninear, int ivar) {
    nopen = nbo[1] + nbo[2];
    const double cnst = 1e-3;
    const int foff = (ivar - 1) * minear, doff = (ivar - 1) * ninear;
    int l = 1, nend = 0, n2 = 0, ninit = 0;
    for (int loop = 1; loop <= 3; ++loop) {
        ninit = nend + 1;
        nend = nend + nbo[loop];
        int n1 = std::max(ninit, nelec + 1);
        n2 = std::min(nend, nelec + nmos);
        if (n2 < n1) continue;
        for (int i = n1; i <= n2; ++i) {
            if (i > ninit) {
                for (int j = ninit; j <= i - 1; ++j) {
                    double diffe = e[i] - e[j];
                    double com = (std::abs(diffe) > 1e-4)
                                     ? (fd[doff + l] - fci[doff + l]) / diffe
                                     : 0.0;
                    cmo[i - 1][j - 1] = -com;
                    cmo[j - 1][i - 1] = com;
                    l++;
                }
            }
            cmo[i - 1][i - 1] = 0.0;
        }
    }
    int ncol = n2 - ninit + 1;
    if (ncol > 0 && n2 < norbs) {
        for (int j = ninit; j <= n2; ++j)
            for (int i = n2 + 1; i <= norbs; ++i) {
                double diffe = e[i] - e[j];
                double com = (std::abs(diffe) > 1e-4)
                                 ? (fd[doff + l] - fci[doff + l]) / diffe
                                 : 0.0;
                cmo[i - 1][j - 1] = -com;
                cmo[j - 1][i - 1] = com;
                l++;
            }
    }
    for (int k = 0; k < nmos; ++k) emo[nelec + 1 + k] = fci[doff + l + k];

    l = 1;
    if (nbo[2] > 0 && nbo[1] > 0) {
        double scal = 1.0 / (2.0 - fract + cnst);
        for (int j = 1; j <= nbo[1]; ++j)
            for (int i = nbo[1] + 1; i <= nopen; ++i) {
                double com = f[foff + l] * scal;
                cmo[i - 1][j - 1] = -com; cmo[j - 1][i - 1] = com; l++;
            }
    }
    if (nbo[3] > 0 && nbo[1] > 0) {
        double scal = 0.5;
        for (int j = 1; j <= nbo[1]; ++j)
            for (int i = nopen + 1; i <= norbs; ++i) {
                double com = f[foff + l] * scal;
                cmo[i - 1][j - 1] = -com; cmo[j - 1][i - 1] = com; l++;
            }
    }
    if (nbo[3] != 0 && nbo[2] != 0) {
        double scal = 1.0 / (fract + cnst);
        for (int j = nbo[1] + 1; j <= nopen; ++j)
            for (int i = nopen + 1; i <= norbs; ++i) {
                double com = f[foff + l] * scal;
                cmo[i - 1][j - 1] = -com; cmo[j - 1][i - 1] = com; l++;
            }
    }
}
