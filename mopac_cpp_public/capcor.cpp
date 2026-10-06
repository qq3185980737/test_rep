// capcor.cpp — C++ translation of MOPAC 2016 "capcor.F90".

#include "capcor.h"

#include "molkst_C.h"

double capcor(const std::vector<int>& nat, const std::vector<int>& nfirst,
              const std::vector<int>& nlast, const std::vector<double>& p,
              const std::vector<double>& h) {
    double sum = 0.0;
    for (int i = 1; i <= molkst_C::numat; ++i) {
        int ni = nat[i];
        int il = nfirst[i];
        int iu = nlast[i];
        if (ni == 102) {
            int j = (nlast[i] * (nlast[i] + 1)) / 2;
            int ii = iu - 1;
            for (int k = 1; k <= ii; ++k) {
                --j;
                sum += p[j] * h[j];
            }
        } else {
            for (int j = 1; j <= i; ++j) {
                int jl = nfirst[j];
                if (nat[j] != 102) continue;
                for (int k = il; k <= iu; ++k) {
                    int kk = (k * (k - 1)) / 2 + jl;
                    sum += p[kk] * h[kk];
                }
            }
        }
    }
    return -2.0 * sum;
}
