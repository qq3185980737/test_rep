// fock2z.cpp — C++ translation of MOPAC 2016 "fock2z.F90".

#include "fock2z.h"

void focd2z(int iab, int jba,
            std::vector<double>& fii, std::vector<double>& fjj,
            std::vector<double>& fij,
            const std::vector<double>& pii, const std::vector<double>& pjj,
            const std::vector<double>& pij,
            const std::vector<double>& wj, const std::vector<double>& wk,
            bool diagonal, int& kr) {
    int loop = 0;
    for (int i = 1; i <= iab; ++i) {
        int ka = (i * (i - 1)) / 2;
        double aa = 2.0;
        for (int j = 1; j <= i; ++j) {
            if (i == j) aa = 1.0;
            int ij = ka + j;
            for (int k = 1; k <= jba; ++k) {
                int kc = (k * (k - 1)) / 2;
                double bb = 2.0;
                for (int l = 1; l <= k; ++l) {
                    if (k == l) bb = 1.0;
                    int kl = kc + l;
                    ++loop;
                    double a = wj[loop];
                    fii[ij] += bb * a * pjj[kl];
                    if (!diagonal) {
                        fjj[kl] += aa * a * pii[ij];
                        a = wk[loop] * aa * bb * 0.125;
                        int ik = (i - 1) * jba + k;
                        int il = (i - 1) * jba + l;
                        int jk = (j - 1) * jba + k;
                        int jl = (j - 1) * jba + l;
                        fij[ik] -= a * pij[jl];
                        fij[il] -= a * pij[jk];
                        fij[jk] -= a * pij[il];
                        fij[jl] -= a * pij[ik];
                    }
                }
            }
        }
    }
    kr += loop;
}

void fock2z(std::vector<double>&, std::vector<double>&,
            std::vector<double>&, const std::vector<double>&,
            const std::vector<double>&,
            std::vector<std::vector<double>>&, int, int) {}
