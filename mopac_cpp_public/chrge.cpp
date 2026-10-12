// chrge.cpp — C++ translation of MOPAC 2016 "chrge.F90".

#include "chrge.h"

#include <vector>

#include "common_arrays_C.h"
#include "ijbo.h"
#include "MOZYME_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace MOZYME_C;

namespace {
    void chrge_for_MOZYME(const std::vector<double>& p, std::vector<double>& q) {
        // F90 chrge_for_MOZYME.F90: q(i) = sum_j p(ijbo(i,i)+j)
        for (int i = 1; i <= numat; ++i) {
            int ii = ijbo(i, i);
            double sum = 0.0;
            for (int j = 1; j <= iorbs[i]; ++j) {
                ii = ii + j;
                sum += p[ii];
            }
            q[i] = sum;
        }
    }
}

void chrge(const std::vector<double>& pin, std::vector<double>& qout) {
    if (mozyme) { chrge_for_MOZYME(pin, qout); return; }
    int k = 0;
    for (int i = 1; i <= numat; ++i) {
        int ia = nfirst[i], ib = nlast[i];
        qout[i] = 0.0;
        for (int j = ia; j <= ib; ++j) {
            k += j;
            qout[i] += pin[k];
        }
    }
}
