// check.cpp — C++ translation of MOPAC 2016 "check.F90".

#include "check.h"

#include <cmath>

#include "MOZYME_C.h"

void check(int nvec, const std::vector<int>& nnc, const std::vector<int>& nc,
           const std::vector<int>& icvec, const std::vector<int>& iorbs,
           const std::vector<int>& ncvec, std::vector<double>& cvec) {
    double error = 0.0;
    for (int i = 1; i <= nvec; ++i) {
        int m = nnc[i];
        double sum = 0.0;
        int n = 0;
        for (int j = 1; j <= nc[i]; ++j) {
            ++m;
            int k = icvec[m];
            for (int mm = 1; mm <= iorbs[k]; ++mm) {
                ++n;
                double v = cvec[ncvec[i] + n];
                sum += v * v;
            }
        }
        error += std::abs(1.0 - sum);
        MOZYME_C::ws[i] = sum;
    }
    for (int i = 1; i <= nvec; ++i) {
        double scale = 1.0 / std::sqrt(MOZYME_C::ws[i]);
        int n = 0;
        for (int j = 1; j <= nc[i]; ++j) {
            int k = icvec[nnc[i] + j];
            for (int mm = 1; mm <= iorbs[k]; ++mm) {
                ++n;
                cvec[ncvec[i] + n] *= scale;
            }
        }
    }
    if (error <= 0.1) return;
    // Severe error: would print and call mopend. In translation, return.
}
