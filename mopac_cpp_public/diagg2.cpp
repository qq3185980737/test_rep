// diagg2.cpp — C++ translation of MOPAC 2016 "diagg2.F90".
//
// DIAGG2 performs a simple Jacobian annihilation of the energy terms
// connecting the occupied LMOS and the virtual LMOS.  The energy terms are
// in the array FMO, and the indices of the LMOS are in IFMO.  With adaptive
// expansion of LMO atom lists (and rejection/retry when array bounds would
// be exceeded).
//
// Signature follows addhb.cpp's forward declaration:
//   eigv:  0-based double* over (nvirtual) eigenvalues (1-based -> data())
//   iused: int* (numat),  latoms: char* (numat, 0/1)
//   storei/storej: double* (norbs) scratch.

#include "diagg2.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace MOZYME_C;
using namespace molkst_C;
using namespace parameters_C;

double reada(const std::string& string, int istart);
namespace {
void epseta(double& eps, double& eta) {
    eps = std::max(std::numeric_limits<double>::min(), 1.0e-39);
    eta = std::max(std::numeric_limits<double>::epsilon(), 1.0e-17);
}
void mopend(const char*) {}
void timer(const char*) {}
}

void diagg2(int nocc, int nvir, const double* eigv, int* iused,
            char* latoms, int nij, int idiagg, double* storei,
            double* storej) {
    static bool debug = false, times = false;
    static int icalcn = 0;
    static double cnst, eps, eta, bigeps;
    static int nrejct[2] = {0, 0};
    if (numcal != icalcn) {
        icalcn = numcal;
        times = (keywrd.find(" TIMES") != std::string::npos);
        debug = (keywrd.find(" DIAGG2") != std::string::npos);
        //
        //   IF THE SYSTEM IS A SOLID, THEN DAMP ROTATION OF VECTORS,
        //   IN AN ATTEMPT TO PREVENT AUTOREGENERATIVE CHARGE OSCILLATION.
        //
        int i = (int)keywrd.find(" DAMP");
        double a;
        if (i > 0) {
            cnst = reada(keywrd, i + 5);
        } else if (id == 3) {
            cnst = 0.5;
        } else {
            a = 1.0;               //  If a transition metal, set DAMP to 0.5
            for (int i2 = 1; i2 <= numat; ++i2)
                if (!main_group[nat[i2]]) a = 0.5;
            cnst = a;
        }
        //
        //   EPS IS THE SMALLEST NUMBER WHICH, WHEN ADDED TO 1.D0, IS NOT
        //   EQUAL TO 1.D0
        //
        epseta(eps, eta);
        //
        //   INCREASE EPS TO ALLOW FOR A LOT OF ROUND-OFF
        //
        bigeps = 50.0 * std::sqrt(eps);
    }
    //
    //  RETRY IS .TRUE. IF THE NUMBER OF REJECTED ANNIHILATIONS IS IN THE
    //  RANGE 1 TO 20 AND THE SAME ON TWO ITERATIONS.
    //
    bool retry = (nrejct[0] == nrejct[1] && nrejct[0] != 0 && nrejct[0] < 20);
    double biglim;
    if (idiagg % 5 == 0 || idiagg <= 5) {
        tiny = -1.0;
        biglim = -1.0;
    } else {
        tiny = 0.01 * tiny;
        biglim = bigeps;
    }
    //
    //   DO A CRUDE 2 BY 2 ROTATION TO "ELIMINATE" SIGNIFICANT ELEMENTS
    //
    for (int i = 1; i <= numat; ++i) { iused[i] = -1; latoms[i] = 0; }
    if (debug) {
        std::fprintf(stdout, "\n");
        std::fprintf(stdout, "            SIZE OF OCCUPIED ARRAYS IN DIAGG2\n\n");
        std::fprintf(stdout, "    LMO    NNCF     NCF   SPACE   NCOCC    SIZE   SPACE\n");
        for (int i = 1; i <= nocc - 1; ++i) {
            int l = ncocc[i];
            for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j)
                l += iorbs[icocc[j]];
            std::fprintf(stdout, "%7d%8d%8d%8d%8d%8d%8d\n", i, nncf[i], ncf[i],
                         nncf[i + 1] - nncf[i] - ncf[i], ncocc[i], l - ncocc[i],
                         ncocc[i + 1] - l);
        }
        int i = nocc;
        int l = ncocc[i];
        for (int j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j)
            l += iorbs[icocc[j]];
        std::fprintf(stdout, "%7d%8d%8d%8d%8d%8d%8d\n", i, nncf[i], ncf[i],
                     icocc_dim - nncf[i] - ncf[i], ncocc[i], l - ncocc[i],
                     cocc_dim - l);
        std::fprintf(stdout, "\n            SIZE OF VIRTUAL ARRAYS IN DIAGG2\n\n");
        std::fprintf(stdout, "    LMO    NNCE     NCE   SPACE   NCVIR    SIZE   SPACE\n");
        for (int i = 1; i <= nvir - 1; ++i) {
            l = ncvir[i];
            for (int j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j)
                l += iorbs[icvir[j]];
            std::fprintf(stdout, "%7d%8d%8d%8d%8d%8d%8d\n", i, nnce[i], nce[i],
                         nnce[i + 1] - nnce[i] - nce[i], ncvir[i], l - ncvir[i],
                         ncvir[i + 1] - l);
        }
        i = nvir;
        l = ncvir[i];
        for (int j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j)
            l += iorbs[icvir[j]];
        std::fprintf(stdout, "%7d%8d%8d%8d%8d%8d%8d\n", i, nnce[i], nce[i],
                     icvir_dim - nnce[i] - nce[i], ncvir[i], l - ncvir[i],
                     cvir_dim - l);
        if (debug && false) {  // 'bug' flag is .false. in source; keep branch dead
            std::fprintf(stdout, "\n");
            std::fprintf(stdout, " THIS FAULT CAN PROBABLY BE CORRECTED BY USE OF KEYWORD 'NLMO=%d'\n\n",
                         ipad2 + 50);
            mopend("VALUE OF NLMO IS TOO SMALL");
        }
    }
    sumb = 0.0;
    int nrej = 0;
    int lij = 0;
    for (int ij = 1; ij <= nij; ++ij) {
        int i = ifmo[1][ij];
        int j = ifmo[2][ij];
        if (std::fabs(fmo[ij]) >= tiny) {
            double c = fmo[ij] * cnst;
            double d = eigs[j] - eigv[i] - shift;
            if (std::fabs(c / d) >= biglim) {
                int ncfj = ncf[j];
                int ncei = nce[i];
                //
                //  STORE LMOS FOR POSSIBLE REJECTION, IF LMOS EXPAND TOO MUCH.
                //
                int jlr = ncocc[j] + 1;
                int jur, jncf;
                if (j != nocc) { jur = ncocc[j + 1]; jncf = nncf[j + 1]; }
                else { jur = cocc_dim; jncf = icocc_dim; }
                jur = std::min(jlr + norbs - 1, jur);
                int ilr = ncvir[i] + 1;
                int iur, incv;
                if (i != nvir) { iur = ncvir[i + 1]; incv = nnce[i + 1]; }
                else { iur = cvir_dim; incv = icvir_dim; }
                iur = std::min(ilr + norbs - 1, iur);
                int l = 0;
                for (int k = jlr; k <= jur; ++k) storej[++l] = cocc[k];
                l = 0;
                for (int k = ilr; k <= iur; ++k) storei[++l] = cvir[k];
                //
                //   STORAGE DONE.
                //
                ++lij;
                double e = std::copysign(std::sqrt(4.0 * c * c + d * d), d);
                double alpha = std::sqrt(0.5 * (1.0 + d / e));
                double beta;
                for (;;) {
                    beta = -std::copysign(std::sqrt(1.0 - alpha * alpha), c);
                    sumb += std::fabs(beta);
                    //
                    // IDENTIFY THE ATOMS IN THE OCCUPIED LMO.  ATOMS NOT USED
                    // ARE FLAGGED BY '-1' IN IUSED.
                    //
                    int mlf = 0;
                    for (int lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf) {
                        int ii = icocc[lf];
                        iused[ii] = mlf;
                        mlf += iorbs[ii];
                    }
                    int loopi = ncvir[i];
                    int loopj = ncocc[j];
                    int mle = 0;
                    //
                    //      ROTATION OF PSEUDO-EIGENVECTORS
                    //
                    for (int le = nnce[i] + 1; le <= nnce[i] + nce[i]; ++le) {
                        int mie = icvir[le];
                        latoms[mie] = 1;
                        int mlff = iused[mie] + loopj;
                        if (iused[mie] >= 0) {
                            //
                            //  TWO BY TWO ROTATION OF ATOMS WHICH ARE COMMON
                            //  TO OCCUPIED LMO J AND VIRTUAL LMO I
                            //
                            for (int mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee) {
                                ++mlff;
                                double a = cocc[mlff];
                                double b = cvir[mlee];
                                cocc[mlff] = alpha * a + beta * b;
                                cvir[mlee] = alpha * b - beta * a;
                            }
                        } else {
                            //
                            //   FILLED LMO ATOM 'MIE' DOES NOT EXIST.
                            //   CHECK IF IT SHOULD EXIST
                            //
                            double sum = 0.0;
                            for (int mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee)
                                sum += (beta * cvir[mlee]) * (beta * cvir[mlee]);
                            if (sum > thresh) {
                                if (nncf[j] + ncf[j] >= jncf) goto reject;
                                if (mlf + iorbs[mie] + loopj > jur) goto reject;
                                //
                                //  YES, OCCUPIED LMO ATOM 'MIE' SHOULD EXIST
                                //
                                ++ncf[j];
                                icocc[nncf[j] + ncf[j]] = mie;
                                iused[mie] = mlf;
                                mlf += iorbs[mie];
                                //
                                //   PUT INTENSITY INTO OCCUPIED LMO ATOM 'MIE'
                                //
                                mlff = iused[mie] + loopj;
                                for (int mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee) {
                                    ++mlff;
                                    cocc[mlff] = beta * cvir[mlee];
                                    cvir[mlee] = alpha * cvir[mlee];
                                }
                            }
                        }
                        mle += iorbs[mie];
                    }
                    //
                    //  NOW CHECK ALL ATOMS WHICH WERE IN THE OCCUPIED LMO
                    //  WHICH ARE NOT IN THE VIRTUAL LMO, TO SEE IF THEY
                    //  SHOULD BE IN THE VIRTUAL LMO.
                    //
                    for (int lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf) {
                        int ii = icocc[lf];
                        if (!latoms[ii]) {
                            double sum = 0.0;
                            for (int mlff = iused[ii] + loopj + 1; mlff <= iused[ii] + loopj + iorbs[ii]; ++mlff)
                                sum += (beta * cocc[mlff]) * (beta * cocc[mlff]);
                            if (sum > thresh) {
                                if (nnce[i] + nce[i] >= incv) goto reject;
                                if (mle + iorbs[ii] + loopi > iur) goto reject;
                                //
                                //  YES, VIRTUAL LMO ATOM 'II' SHOULD EXIST
                                //
                                ++nce[i];
                                icvir[nnce[i] + nce[i]] = ii;
                                latoms[ii] = 1;
                                //
                                //   PUT INTENSITY INTO VIRTUAL LMO ATOM 'II'
                                //
                                int mlff = iused[ii] + loopj;
                                for (int mlee = mle + 1 + loopi; mlee <= mle + iorbs[ii] + loopi; ++mlee) {
                                    ++mlff;
                                    cvir[mlee] = -beta * cocc[mlff];
                                    cocc[mlff] = alpha * cocc[mlff];
                                }
                                mle += iorbs[ii];
                            }
                        }
                    }
                    break;
                reject:
                    ++nrej;
                    //
                    //   THE ARRAY BOUNDS WERE GOING TO BE EXCEEDED.
                    //   TO PREVENT THIS, RESET THE LMOS.
                    //
                    l = 0;
                    for (int k = jlr; k <= jur; ++k) cocc[k] = storej[++l];
                    l = 0;
                    for (int k = ilr; k <= iur; ++k) cvir[k] = storei[++l];
                    ncf[j] = ncfj;
                    nce[i] = ncei;
                    for (int k = 1; k <= numat; ++k) { iused[k] = -1; latoms[k] = 0; }
                    if (retry) {
                        //
                        //   HALF THE ROTATION ANGLE.  WILL THIS PREVENT THE
                        //   ARRAY BOUND FROM BEING EXCEEDED?
                        //
                        alpha = 0.5 * (alpha + 1.0);
                    } else {
                        goto next_ij;
                    }
                }
                //
                //  RESET COUNTERS WHICH HAVE BEEN SET.
                //
                for (int le = nnce[i] + 1; le <= nnce[i] + nce[i]; ++le)
                    latoms[icvir[le]] = 0;
                for (int lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf)
                    iused[icocc[lf]] = -1;
            }
        }
    next_ij: ;
    }
    nrejct[1] = nrejct[0];
    nrejct[0] = nrej;
    if (times) timer(" AFTER DIAGG2 IN ITER");
}
