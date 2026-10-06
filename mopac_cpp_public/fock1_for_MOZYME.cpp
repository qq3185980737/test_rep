// fock1_for_MOZYME.cpp

#include "fock1_for_MOZYME.h"

#include <algorithm>

void fock1_for_MOZYME(std::vector<double>& f, const std::vector<double>& ptot,
                      const std::vector<std::vector<double>>& w,
                      int& kr, int iab, int ilim) {
    for (int i = 1; i <= iab; ++i)
        for (int j = 1; j <= i; ++j) {
            int ij = (i * (i - 1)) / 2 + j;
            double sum = 0.0;
            for (int k = 1; k <= iab; ++k)
                for (int l = 1; l <= iab; ++l) {
                    int ip = std::max(k, l), jp = std::min(k, l);
                    int ijp = (ip * (ip - 1)) / 2 + jp;
                    int klw = ijp;
                    int im = std::max(k, j), jm = std::min(k, j);
                    int ikw = (im * (im - 1)) / 2 + jm;
                    im = std::max(l, i); jm = std::min(l, i);
                    int jlw = (im * (im - 1)) / 2 + jm;
                    sum += ptot[ijp] * w[ij][klw] - 0.5 * ptot[ijp] * w[ikw][jlw];
                }
            f[ij] += sum;
        }
    kr += ilim * ilim;
}
