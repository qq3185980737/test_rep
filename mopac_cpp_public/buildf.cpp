// buildf.cpp — C++ translation of MOPAC 2016 "buildf.F90".

#include "buildf.h"

#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

namespace {
void memory_error(const char*) {}
// fock2z: two-electron Fock construction (not ported). Stub.
void fock2z(std::vector<double>&, std::vector<double>&, std::vector<double>&,
            std::vector<double>&, std::vector<double>&,
            std::vector<std::vector<double>>&, int, int) {}
}

void buildf(std::vector<double>& F, const std::vector<double>& partf, int mode) {
    std::vector<double> q(numat + 1, 0.0), qe(numat + 1, 0.0);
    std::vector<std::vector<double>> ptot2(
        numat + 1, std::vector<double>(82, 0.0));

    if (mode == -1) {
        for (int k = 1; k <= mpack; ++k) F[k] = partf[k] - h[k];
    } else if (mode == 0) {
        for (int k = 1; k <= mpack; ++k) F[k] = h[k];
    } else {  // mode == 1
        for (int k = 1; k <= mpack; ++k) F[k] = partf[k] + h[k];
    }
    if (id == 0) {
        fock2z(F, q, qe, w, w, ptot2, mode, 1);
    } else {
        fock2z(F, q, qe, w, wk, ptot2, mode, 0);
    }
}
