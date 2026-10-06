// adjvec.cpp — C++ translation of MOPAC 2016 "adjvec.F90" (103 lines).
//
// Fortran 1-based indexing preserved: all packed vectors are accessed with
// their original 1-based subscripts (vector index 0 is padding).

#include "adjvec.h"

#include <cmath>

#include "MOZYME_C.h"
#include "molkst_C.h"

void adjvec(std::vector<double>& cvecb, int ncvb,
            std::vector<int>& icvecb, int nib,
            const std::vector<int>& nncb, std::vector<int>& ncb_loc, int nnb,
            const std::vector<int>& ncvecb, int lmob,
            const std::vector<int>& iorbs,
            const std::vector<double>& cveca, int /*ncva*/,
            const std::vector<int>& icveca, int /*nia*/,
            const std::vector<int>& nnca, const std::vector<int>& nca_loc,
            int /*nna_loc*/, const std::vector<int>& ncveca, int lmoa,
            double beta, std::vector<int>& iused, double& sumtot) {
    int ii=0, la=0, lb=0, llim=0, mia=0, mla=0, mlb=0, mlim=0, mlla=0, mllb=0;
    double cutoff, sum;

    cutoff = MOZYME_C::thresh * 1.0e1;
    if (std::fabs(beta) < cutoff) return;
    sumtot = sumtot + std::fabs(beta);

    // Flag atoms of LMO-a that are not (yet) in LMO-b with -1.
    for (la = nnca[lmoa] + 1; la <= nnca[lmoa] + nca_loc[lmoa]; ++la) {
        iused[icveca[la]] = -1;
    }
    mlb = ncvecb[lmob];
    if (lmob == nnb) llim = nib;
    else llim = nncb[lmob + 1];

    if (lmob == nnb) mlim = ncvb - 4;
    else mlim = ncvecb[lmob + 1] - 4;

    for (lb = nncb[lmob] + 1; lb <= nncb[lmob] + ncb_loc[lmob]; ++lb) {
        ii = icvecb[lb];
        iused[ii] = mlb;
        mlb = mlb + iorbs[ii];
    }
    mla = ncveca[lmoa];

    // Rotate the second vector to make it orthogonal to the first.
    for (la = nnca[lmoa] + 1; la <= nnca[lmoa] + nca_loc[lmoa]; ++la) {
        mia = icveca[la];
        if (iused[mia] >= 0) {
            mllb = iused[mia];
            // atoms common to both LMOs
            for (mlla = mla + 1; mlla <= mla + iorbs[mia]; ++mlla) {
                mllb = mllb + 1;
                cvecb[mllb] = cvecb[mllb] - beta * cveca[mlla];
            }
        } else {
            // atom 'mia' not in lmob; if it should be, make it exist.
            sum = 0.0;
            for (mlla = mla + 1; mlla <= mla + iorbs[mia]; ++mlla) {
                sum = sum + cveca[mlla] * cveca[mlla];
            }
            if (beta * beta * sum > cutoff) {
                if (ncb_loc[lmob] < llim && mlb < mlim) {
                    ncb_loc[lmob] = ncb_loc[lmob] + 1;
                    icvecb[nncb[lmob] + ncb_loc[lmob]] = mia;
                    iused[mia] = mlb;
                    for (mlla = mla + 1; mlla <= mla + iorbs[mia]; ++mlla) {
                        mlb = mlb + 1;
                        cvecb[mlb] = -beta * cveca[mlla];
                    }
                }
            }
        }
        mla = mla + iorbs[mia];
    }
    if (mlla != -1) return;

    // Debug block: sum of overlap after orthogonalisation (normally skipped).
    sum = 0.0;
    mla = ncveca[lmoa];
    for (la = nnca[lmoa] + 1; la <= nnca[lmoa] + nca_loc[lmoa]; ++la) {
        mia = icveca[la];
        if (iused[mia] >= 0) {
            mllb = iused[mia];
            for (mlla = mla + 1; mlla <= mla + iorbs[mia]; ++mlla) {
                mllb = mllb + 1;
                sum = sum + cvecb[mllb] * cveca[mlla];
            }
        }
        mla = mla + iorbs[mia];
    }
}
