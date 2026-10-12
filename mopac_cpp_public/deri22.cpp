// deri22.cpp — C++ translation of MOPAC 2016 "deri22.F90" (complete).
// Builds the supervector AB = (DIAG + C'*FOC2*C)*B used by the relaxation
// step of deri2, plus the CI-active Fock diagonal blocks FCI.
//
// Layout notes: c is [row][col] 1-based; internal flat cf/wflat/dp use the
// 1-based-padding convention (physical offset == Fortran element number,
// element 0 unused). b/ab/fci are column-major supervectors over the full
// basis set; bcol/abcol/fcicol select the current 1-based column so the
// caller passes the same container for every column (as F90 b(1,j)).
#include "deri22.h"

#include <algorithm>
#include <vector>

#include "common_arrays_C.h"
#include "fock2.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mxm.h"
#include "mxmt.h"
#include "mtxm.h"
#include "supdot.h"

using namespace common_arrays_C;
using namespace meci_C;
using namespace molkst_C;

void deri22(const std::vector<std::vector<double>>& c,
            std::vector<double>& b, int bcol,
            std::vector<std::vector<double>>& work,
            std::vector<double>& foc2,
            std::vector<double>& ab, int abcol, int minear,
            std::vector<double>& fci, int fcicol,
            std::vector<double>& w,
            const std::vector<double>& diag,
            const std::vector<double>& scalar, int ninear) {
    (void)ninear;
    std::vector<double> cf(norbs * norbs + 1, 0.0);
    for (int i = 1; i <= norbs; ++i)
        for (int ip = 1; ip <= norbs; ++ip)
            cf[(i - 1) * norbs + ip] = c[ip][i];
    std::vector<double> dp(norbs * norbs + 1, 0.0);
    std::vector<double> dpa(norbs * norbs + 1, 0.0);
    std::vector<double> wflat(norbs * norbs + 1, 0.0);
    const int boff = (bcol - 1) * minear, aoff = (abcol - 1) * minear,
              foff = (fcicol - 1) * ninear;

    // STEP 0: unscale B.
    for (int i = 1; i <= minear; ++i) b[boff + i] *= scalar[i];

    // STEP 1: work = C*B in the packed block order of B.
    int l = 1;
    if (nbo[2] != 0 && nbo[1] != 0) {
        // OPEN-CLOSED
        mxm(&cf[nbo[1] * norbs + 1], norbs, &b[boff + l], nbo[2],
            &wflat[1], nbo[1]);
        // CLOSED-OPEN
        mxmt(&cf[1], norbs, &b[boff + l], nbo[1],
             &wflat[nbo[1] * norbs + 1], nbo[2]);
        l += nbo[2] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[1] != 0) {
        // VIRTUAL-CLOSED
        if (l > 1) {
            mxm(&cf[nopen * norbs + 1], norbs, &b[boff + l], nbo[3],
                &dp[1], nbo[1]);
            int icount = 0;
            for (int j = 1; j <= nbo[1]; ++j)
                for (int i = 1; i <= norbs; ++i) {
                    ++icount;
                    wflat[(j - 1) * norbs + i] += dp[icount];
                }
        } else {
            mxm(&cf[nopen * norbs + 1], norbs, &b[boff + l], nbo[3],
                &wflat[1], nbo[1]);
        }
        // CLOSED-VIRTUAL
        mxmt(&cf[1], norbs, &b[boff + l], nbo[1],
             &wflat[nopen * norbs + 1], nbo[3]);
        l += nbo[3] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[2] != 0) {
        // VIRTUAL-OPEN
        mxm(&cf[nopen * norbs + 1], norbs, &b[boff + l], nbo[3],
            &dp[1], nbo[2]);
        {
            int icount = 0;
            for (int j = nbo[1] + 1; j <= nbo[1] + nbo[2]; ++j)
                for (int i = 1; i <= norbs; ++i) {
                    ++icount;
                    wflat[(j - 1) * norbs + i] += dp[icount];
                }
        }
        // OPEN-VIRTUAL
        mxmt(&cf[nbo[1] * norbs + 1], norbs, &b[boff + l], nbo[2],
             &dp[1], nbo[3]);
        {
            int icount = 0;
            for (int j = nopen + 1; j <= nopen + nbo[3]; ++j)
                for (int i = 1; i <= norbs; ++i) {
                    ++icount;
                    wflat[(j - 1) * norbs + i] += dp[icount];
                }
        }
    }


    // STEP 2: dp = work*C', packed canonical.
    int l2 = 0;
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= i; ++j) {
            double s = 0.0;
            for (int m = 1; m <= norbs; ++m)
                s += wflat[(m - 1) * norbs + i] * cf[(m - 1) * norbs + j];
            dp[++l2] = s;
        }

    // 2-electron Fock matrix built with the density-matrix derivative.
    std::fill(foc2.begin(), foc2.end(), 0.0);
    for (int i = 1; i <= norbs * norbs; ++i) dpa[i] = 0.5 * dp[i];
    std::vector<double> wj(2, 0.0), wk(2, 0.0);
    fock2(foc2, dp, dpa, w, wj, wk, numat, nfirst, nlast, 2);
    // dp(norbs,nend) = foc2(norbs,norbs) * c(norbs,nend).
    int nend = std::max(nopen, nelec + nmos);
    for (int i = 1; i <= nend; ++i)
        supdot(&dp[(i - 1) * norbs], foc2.data(), &cf[(i - 1) * norbs],
               norbs, 1);

    // EXTRACT FCI.
    l = 1;
    int nend2 = 0, ninit = 0, n2 = 0;
    for (int loop = 1; loop <= 3; ++loop) {
        ninit = nend2 + 1;
        nend2 = nend2 + nbo[loop];
        int n1 = std::max(ninit, nelec + 1);
        n2 = std::min(nend2, nelec + nmos);
        if (n2 < n1) continue;
        for (int i = n1; i <= n2; ++i) {
            if (i <= ninit) continue;
            mxm(&cf[(i - 1) * norbs + 1], 1,
                &dp[(ninit - 1) * norbs + 1], norbs, &fci[foff + l],
                i - ninit);
            l += i - ninit;
        }
    }
    int ncol = n2 - ninit + 1;
    if (ncol > 0 && n2 < norbs) {
        mtxm(&cf[n2 * norbs + 1], norbs - n2,
             &dp[(ninit - 1) * norbs + 1], norbs, &fci[foff + l], ncol);
        l += ncol * (norbs - n2);
    }
    for (int i = nelec + 1; i <= nelec + nmos; ++i) {
        double s = 0.0;
        for (int k = 1; k <= norbs; ++k)
            s += cf[(i - 1) * norbs + k] * dp[(i - 1) * norbs + k];
        fci[foff + l] = -s;
        ++l;
    }

    // NEW SUPERVECTOR AB = (DIAG + C'*FOC2*C)*B, SCALED.
    l = 1;
    if (nbo[2] != 0 && nbo[1] != 0) {
        mtxm(&cf[nbo[1] * norbs + 1], nbo[2], &dp[1], norbs,
             &ab[aoff + l], nbo[1]);
        l += nbo[2] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[1] != 0) {
        mtxm(&cf[nopen * norbs + 1], nbo[3], &dp[1], norbs,
             &ab[aoff + l], nbo[1]);
        l += nbo[3] * nbo[1];
    }
    if (nbo[3] != 0 && nbo[2] != 0)
        mtxm(&cf[nopen * norbs + 1], nbo[3],
             &dp[norbs * nbo[1] + 1], norbs, &ab[aoff + l], nbo[2]);

    // PART 2: ab = scale*(diag*b + ab); rescale basis vector b.
    for (int i = 1; i <= minear; ++i) {
        ab[aoff + i] = (diag[i] * b[boff + i] + ab[aoff + i]) * scalar[i];
        b[boff + i] = b[boff + i] / scalar[i];
    }

    // Copy the flat work buffer back into the out-parameter (1-based rows).
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= norbs; ++j)
            work[i][j] = wflat[(j - 1) * norbs + i];
}
