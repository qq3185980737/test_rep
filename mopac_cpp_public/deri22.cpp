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
    fprintf(stderr, "[D22BS] %+15.8f %+15.8f %+15.8f %+15.8f %+15.8f %+15.8f %+15.8f %+15.8f\n",
        b[boff+1],b[boff+2],b[boff+3],b[boff+4],b[boff+5],b[boff+6],b[boff+7],b[boff+8]); fflush(stderr);
    if (minear > 8)
        fprintf(stderr, "[D22B9] %+15.8f %+15.8f %+15.8f %+15.8f\n",
            b[boff+9],b[boff+10],b[boff+11],b[boff+12]); fflush(stderr);
    fprintf(stderr, "[D22C5] %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f\n",
        cf[1],cf[2],cf[3],cf[4],cf[5],
        cf[norbs+1],cf[norbs+2],cf[norbs+3],cf[norbs+4],cf[norbs+5],
        cf[2*norbs+1],cf[2*norbs+2],cf[2*norbs+3],cf[2*norbs+4],cf[2*norbs+5],
        cf[3*norbs+1],cf[3*norbs+2],cf[3*norbs+3],cf[3*norbs+4],cf[3*norbs+5],
        cf[4*norbs+1],cf[4*norbs+2],cf[4*norbs+3],cf[4*norbs+4],cf[4*norbs+5]); fflush(stderr);

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

    fprintf(stderr, "[D22W] %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
        wflat[1],wflat[2],wflat[3],wflat[4],wflat[5],
        wflat[norbs+1],wflat[norbs+2],wflat[norbs+3],wflat[norbs+4],wflat[norbs+5],
        wflat[2*norbs+1],wflat[2*norbs+2],wflat[2*norbs+3],wflat[2*norbs+4],wflat[2*norbs+5],
        wflat[3*norbs+1],wflat[3*norbs+2],wflat[3*norbs+3],wflat[3*norbs+4],wflat[3*norbs+5],
        wflat[4*norbs+1],wflat[4*norbs+2],wflat[4*norbs+3],wflat[4*norbs+4],wflat[4*norbs+5]); fflush(stderr);

    // STEP 2: dp = work*C', packed canonical.
    int l2 = 0;
    for (int i = 1; i <= norbs; ++i)
        for (int j = 1; j <= i; ++j) {
            double s = 0.0;
            for (int m = 1; m <= norbs; ++m)
                s += wflat[(m - 1) * norbs + i] * cf[(m - 1) * norbs + j];
            dp[++l2] = s;
        }
    fprintf(stderr, "[D22DP] %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f %+12.7f\n",
        dp[1],dp[2],dp[3],dp[4],dp[5],dp[6],dp[7],dp[8],dp[9],dp[10],dp[11],dp[12],dp[13],dp[14],dp[15]); fflush(stderr);

    // 2-electron Fock matrix built with the density-matrix derivative.
    std::fill(foc2.begin(), foc2.end(), 0.0);
    for (int i = 1; i <= norbs * norbs; ++i) dpa[i] = 0.5 * dp[i];
    std::vector<double> wj(2, 0.0), wk(2, 0.0);
    fock2(foc2, dp, dpa, w, wj, wk, numat, nfirst, nlast, 2);
    fprintf(stderr, "[D22F2] %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n",
        foc2[1],foc2[2],foc2[3],foc2[4],foc2[5],foc2[6],foc2[7],foc2[8],foc2[9],foc2[10],
        foc2[11],foc2[12],foc2[13],foc2[14],foc2[15],foc2[16],foc2[17],foc2[18],foc2[19],foc2[20],foc2[21]); fflush(stderr);
    // dp(norbs,nend) = foc2(norbs,norbs) * c(norbs,nend).
    int nend = std::max(nopen, nelec + nmos);
    for (int i = 1; i <= nend; ++i)
        supdot(&dp[(i - 1) * norbs], foc2.data(), &cf[(i - 1) * norbs],
               norbs, 1);
    fprintf(stderr, "[D22DN] %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f %+.6f\n",
        dp[1],dp[2],dp[3],dp[4],dp[5],dp[6],dp[7],dp[8],dp[9],dp[10],dp[11],dp[12],dp[13],dp[14],dp[15],dp[16],dp[17],dp[18],dp[19],dp[20],dp[21],dp[22],dp[23],dp[24],dp[25],dp[26],dp[27],dp[28],dp[29],dp[30],dp[31],dp[32],dp[33],dp[34],dp[35],dp[36]); fflush(stderr);

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
