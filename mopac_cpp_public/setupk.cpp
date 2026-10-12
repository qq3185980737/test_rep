// setupk.cpp — C++ translation of MOPAC 2016 "setupk.F90".
// Uses MOZYME_C (shared) arrays icocc/ncf/nncf/kopt (1-based).
#include "setupk.h"
#include <algorithm>
#include "molkst_C.h"
#include "MOZYME_C.h"
using molkst_C::numat;
using namespace MOZYME_C;
void setupk(int nocc1) {
    std::fill(kopt.begin(), kopt.end(), 0);
    for (int i = 1; i <= nocc1; ++i) {
        int j = nncf[i];
        for (int k = 1; k <= ncf[i]; ++k)
            kopt[icocc[k + j]] = 1;
    }
    int l = 0;
    for (int i = 1; i <= numat; ++i) {
        if (kopt[i] == 1) { ++l; kopt[l] = i; }
    }
    if (l != numat) kopt[l + 1] = 0;
}
