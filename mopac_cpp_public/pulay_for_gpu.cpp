// pulay_for_gpu.cpp — Pulay DIIS convergence accelerator (non-GPU path).
#include "pulay_for_gpu.h"

#include <cmath>
#include <vector>

#include "blas1.h"
#include "molkst_C.h"
#include "mult_symm_AB.h"
#include "osinv.h"

using namespace molkst_C;

void pulay_for_gpu(double* f, double* p, int n, double* fppf, double* fock,
                   double* emat, int& lfock, int& nfock, int msize,
                   bool& start, double& pl) {
    static int icalcn = 0;
    static int maxlim = 6;
    static bool debug = false;
    int linear = n * (n + 1) / 2;
    int mfock = msize / linear;
    if (mfock > maxlim) mfock = maxlim;

    if (start) {
        nfock = 1;
        lfock = 1;
        start = false;
    } else {
        if (nfock < mfock) nfock = nfock + 1;
        if (lfock != mfock) lfock = lfock + 1;
        else lfock = 1;
    }
    int lbase = (lfock - 1) * linear;

    // Store F (column-stride mfock).
    for (int k = 0; k < linear; ++k)
        fock[lfock - 1 + k * mfock] = f[k];

    int iopc = lgpu ? 4 : 3;

    // FPPF = P*F - F*P.
    mult_symm_AB(p, f, 1.0, n, linear, fppf + lbase, 0.0, iopc);
    mult_symm_AB(f, p, 1.0, n, linear, fppf + lbase, -1.0, iopc);

    int nfock1 = nfock + 1;
    for (int i = 1; i <= nfock; ++i) {
        emat[(nfock1 - 1) * 20 + i - 1] = -1.0;
        emat[(i - 1) * 20 + nfock1 - 1] = -1.0;
        double dot = 0.0;
        for (int dd = 0; dd < linear; ++dd)
            dot += fppf[(i - 1) * linear + dd] * fppf[lbase + dd];
        emat[(lfock - 1) * 20 + i - 1] = dot;
        emat[(i - 1) * 20 + lfock - 1] = dot;
    }
    pl = emat[(lfock - 1) * 20 + lfock - 1] / linear;
    emat[(nfock1 - 1) * 20 + nfock1 - 1] = 0.0;
    double diag = emat[(lfock - 1) * 20 + lfock - 1];
    if (diag < 1e-20) return;
    double cconst = 1.0 / diag;
    for (int i = 0; i < nfock; ++i)
        for (int j = 0; j < nfock; ++j)
            emat[i * 20 + j] *= cconst;

    // Copy emat into evec (nfock1 x nfock1 row-major).
    std::vector<double> evec(nfock1 * nfock1, 0.0);
    for (int i = 0; i < nfock1; ++i)
        for (int j = 0; j < nfock1; ++j)
            evec[i * nfock1 + j] = emat[i * 20 + j];

    double d = 0.0;
    osinv(evec.data(), nfock1, d);
    if (std::fabs(d) < 1e-6) { start = true; return; }
    if (nfock < 2) return;

    // coeffs = last column of inverse, negated.
    std::vector<double> coeffs(nfock, 0.0);
    for (int i = 0; i < nfock; ++i)
        coeffs[i] = -evec[i * nfock1 + nfock1 - 1];

    // Best Fock = linear combination.
    for (int i = 0; i < linear; ++i) {
        double sum = 0.0;
        int ii = i * mfock;
        for (int j = 0; j < nfock; ++j)
            sum += coeffs[j] * fock[j + ii];
        f[i] = sum;
    }
}
