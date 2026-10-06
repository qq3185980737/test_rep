// chrge.cpp — C++ translation of MOPAC 2016 "chrge.F90".

#include "chrge.h"

#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

namespace { void chrge_for_MOZYME(const std::vector<double>&, std::vector<double>&) {} }

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
