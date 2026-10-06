// molval.cpp
#include "molval.h"
#include <algorithm>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

void molval(double* c, double* p, double rhfuhf) {
    std::vector<double> val(norbs+1, 0.0);
    for (int i = 1; i <= norbs; ++i) {
        double sum = 0.0;
        for (int jj = 1; jj <= numat; ++jj) {
            int jl = nfirst[jj], ju = nlast[jj];
            for (int j = jl; j <= ju; ++j) {
                for (int kk = 1; kk <= numat; ++kk) {
                    if (kk == jj) continue;
                    int kl = nfirst[kk], ku = nlast[kk];
                    for (int k = kl; k <= ku; ++k) {
                        int l1 = std::max(j, k);
                        int l2 = j + k - l1;
                        int l = (l1*(l1-1))/2 + l2;
                        sum += c[(i-1)*norbs + j-1] * c[(i-1)*norbs + k-1] * p[l];
                    }
                }
            }
        }
        val[i] = sum * rhfuhf;
    }
}
