// deri21.cpp — C++ translation of MOPAC 2016 "deri21.F90".
// Least-squares orthonormal basis: eigendecomposition of A'A (packed), then
// B = A * V * diag(|lambda|^-1/2) truncated at cumulative variance cutoff.
// Layout: a/b column-major 0-based (F90 a(minear,nvar_nvo), b(minear,ncut)).
#include "deri21.h"
#include <cmath>
#include <vector>
#include "mtxmc.h"
#include "rsp.h"
#include "mxm.h"

void deri21(const double* a, int nvar_nvo, int minear, int ifirst,
            double* vnert, double* pnert, double* b, int& ncut) {
    const double cutoff = 0.85;
    double sum2 = 0.0;
    // VNERT = A' * A, packed lower triangle (size nvar_nvo*(nvar_nvo+1)/2)
    std::vector<double> work((nvar_nvo * nvar_nvo + 1) * 2, 0.0);
    mtxmc(a, nvar_nvo, a, minear, work.data());
    int np = (nvar_nvo * (nvar_nvo + 1)) / 2;
    for (int k = 0; k < np; ++k) work[k] = -work[k];
    if (std::fabs(work[0]) < 1e-28 && nvar_nvo == 1) {
        pnert[0] = std::sqrt(-work[0]);
        work[0] = 1e15;
        vnert[0] = 1.0;
        ncut = 1;
    } else {
        rsp(work.data(), nvar_nvo, pnert, vnert);
        // eigenvalues ascending (LAPACK convention, matching MOPAC rsp)
        double sum = 0.0;
        for (int i = 0; i < nvar_nvo; ++i) sum -= pnert[i];
        int l = 0;
        int i = 0;
        for (; i < ifirst; ++i) {
            sum2 -= pnert[i] / sum;
            pnert[i] = std::sqrt(std::fabs(pnert[i]));
            for (int q = 0; q < nvar_nvo; ++q) work[l + q] = vnert[l + q] / pnert[i];
            l += nvar_nvo;
            if (sum2 < cutoff) continue;
            ncut = i + 1;
            break;
        }
        if (i >= ifirst) ncut = ifirst;
    }
    // B(minear,ncut) = A(minear,nvar_nvo) * work(nvar_nvo,ncut)
    mxm(a, minear, work.data(), nvar_nvo, b, ncut);
}
