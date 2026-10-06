// reorth.cpp — C++ translation of MOPAC 2016 "reorth.F90".
// Re-orthogonalise the MOZYME LMOs: builds ws (coefficients in atom order)
// for each LMO in turn, computes overlap with every later LMO of the same
// space and every LMO of the other space, and rotates via adjvec. All
// MOZYME bookkeeping arrays are Fortran 1-based (index 0 unused).
#include "reorth.h"

#include <cstdio>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
using namespace MOZYME_C;

using common_arrays_C::nfirst;
// Bare globals (declared outside namespace MOZYME_C).
using MOZYME_C::nvirtual;
using MOZYME_C::noccupied;

void adjvec(std::vector<double>& cvecb, int ncvb,
            std::vector<int>& icvecb, int nib,
            const std::vector<int>& nncb, std::vector<int>& ncb_loc, int nnb,
            const std::vector<int>& ncvecb, int lmob,
            const std::vector<int>& iorbs,
            const std::vector<double>& cveca, int ncva,
            const std::vector<int>& icveca, int nia,
            const std::vector<int>& nnca, const std::vector<int>& nca_loc,
            int nna_loc, const std::vector<int>& ncveca, int lmoa,
            double beta, std::vector<int>& iused, double& sumtot);

void reorth(double* ws) {
    using MOZYME_C::iorbs;
    const int numat = molkst_C::numat;
    std::vector<char> latom_loc(numat + 1, 0);
    std::vector<int> iused(numat + 1, 0);
    double sumtot = 0.0;

    for (int i = 1; i <= nvirtual; ++i) {
        int loopi = ncvir[i];
        for (int j = 1; j <= numat; ++j) latom_loc[j] = 0;
        for (int jj = nnce[i] + 1; jj <= nnce[i] + nce[i]; ++jj) {
            const int j1 = icvir[jj];
            latom_loc[j1] = 1;
            int j = nfirst[j1] - 1;
            for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                ++loopi;
                ++j;
                ws[j] = cvir[loopi];
            }
        }
        for (int ii = i + 1; ii <= nvirtual; ++ii) {
            double sum = 0.0;
            int loopii = ncvir[ii];
            for (int jj = nnce[ii] + 1; jj <= nnce[ii] + nce[ii]; ++jj) {
                const int j1 = icvir[jj];
                if (latom_loc[j1]) {
                    int j = nfirst[j1] - 1;
                    for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                        ++loopii;
                        ++j;
                        sum += ws[j] * cvir[loopii];
                    }
                } else {
                    loopii += iorbs[j1];
                }
            }
            adjvec(cvir, cvir_dim, icvir, icvir_dim, nnce, nce, nvirtual,
                   ncvir, i, iorbs, cvir, cvir_dim, icvir, icvir_dim, nnce,
                   nce, nvirtual, ncvir, ii, sum, iused, sumtot);
        }
        for (int ii = 1; ii <= noccupied; ++ii) {
            double sum = 0.0;
            int loopii = ncocc[ii];
            for (int jj = nncf[ii] + 1; jj <= nncf[ii] + ncf[ii]; ++jj) {
                const int j1 = icocc[jj];
                if (latom_loc[j1]) {
                    int j = nfirst[j1] - 1;
                    for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                        ++loopii;
                        ++j;
                        sum += ws[j] * cocc[loopii];
                    }
                } else {
                    loopii += iorbs[j1];
                }
            }
            adjvec(cvir, cvir_dim, icvir, icvir_dim, nnce, nce, nvirtual,
                   ncvir, i, iorbs, cocc, cocc_dim, icocc, icocc_dim, nncf,
                   ncf, noccupied, ncocc, ii, sum, iused, sumtot);
        }
    }

    for (int i = 1; i <= noccupied; ++i) {
        int loopi = ncocc[i];
        for (int j = 1; j <= numat; ++j) latom_loc[j] = 0;
        for (int jj = nncf[i] + 1; jj <= nncf[i] + ncf[i]; ++jj) {
            const int j1 = icocc[jj];
            latom_loc[j1] = 1;
            int j = nfirst[j1] - 1;
            for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                ++loopi;
                ++j;
                ws[j] = cocc[loopi];
            }
        }
        for (int ii = i + 1; ii <= noccupied; ++ii) {
            double sum = 0.0;
            int loopii = ncocc[ii];
            for (int jj = nncf[ii] + 1; jj <= nncf[ii] + ncf[ii]; ++jj) {
                const int j1 = icocc[jj];
                if (latom_loc[j1]) {
                    int j = nfirst[j1] - 1;
                    for (int jx = 1; jx <= iorbs[j1]; ++jx) {
                        ++loopii;
                        ++j;
                        sum += ws[j] * cocc[loopii];
                    }
                } else {
                    loopii += iorbs[j1];
                }
            }
            adjvec(cocc, cocc_dim, icocc, icocc_dim, nncf, ncf, noccupied,
                   ncocc, i, iorbs, cocc, cocc_dim, icocc, icocc_dim, nncf,
                   ncf, noccupied, ncocc, ii, sum, iused, sumtot);
        }
    }
    std::printf(" Total error:%19.12f\n", sumtot);
}
