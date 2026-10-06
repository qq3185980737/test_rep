// diagg.cpp — C++ translation of diagg.F90 / diagg1.F90 / diagg2.F90 /
// density_for_MOZYME.F90 / epseta.F90 (MOPAC 2016).
#include "diagg.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "ijbo.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace molkst_C;

double reada(const std::string&, int);  // external (1-based position)
void timer(const std::string& text);  // external timing hook

namespace {

// diagg1 statics (Fortran save).
int nf = 0, icalcn1 = 0, mydisp = 0, ij0 = 0;
bool times1 = false;
double fref = 10.0, oldlim = 0.0, safety = 1.0;

// diagg2 statics (Fortran save).
int icalcn2 = 0;
bool debug2 = false, times2 = false;
double const2 = 1.0, eps2 = 1.0, eta2 = 1.0, bigeps = 1.0;
int nrejct[2] = {0, 0};

}  // namespace

void epseta(double& eps, double& eta) {
    // RETURN ETA, THE SMALLEST REPRESENTABLE NUMBER,
    // AND EPS, THE SMALLEST NUMBER FOR WHICH 1+EPS.NE.1.
    eps = std::max(std::numeric_limits<double>::min(), 1.e-39);
    eta = std::max(std::numeric_limits<double>::epsilon(), 1.e-17);
}

void diagg1(const std::vector<double>& fao, int nocc, int nvir,
            std::vector<double>& eigv, std::vector<double>& ws,
            std::vector<bool>& latoms, std::vector<std::vector<int>>& ifmo,
            std::vector<double>& fmo, int fmo_dim, int& nij, int idiagg,
            std::vector<double>& avir, std::vector<double>& aocc,
            std::vector<double>& aov) {
    // eigv is the eigs(nocc+1:) slice: element i (1..nvir) lives at
    // eigv[nocc+i] in the full vector (which is what the caller passes).
    int i, i1, i2, i4, ii, j, j1, j2, j4, jj, jl, jx, k, k1, kk, kl, l, loopi,
        loopj, kj, i1j1, i1j2, i2j1, i2j2, k1j1;
    bool lij;
    double cutlim, flim, cutoff, sum, sum1;
    if (numcal != icalcn1) {
        icalcn1 = numcal;
        fref = 10.0;
        safety = 1.0;
        oldlim = 0.0;
        nf = 0;
        times1 = (keywrd.find(" TIMES") != std::string::npos);
        if (keywrd.find(" OLDENS") != std::string::npos) {
            fref = 0.0;
        }
    }
    std::fill(aocc.begin(), aocc.end(), 0.0);
    //  contributions of atoms in the occupied LMO's -> aocc (square)
    for (j = 1; j <= nocc; ++j) {
        loopj = ncocc[j];
        kl = 0;
        for (kk = nncf[j] + 1; kk <= nncf[j] + ncf[j]; ++kk) {
            k1 = icocc[kk];
            sum = 0.0;
            for (k = nfirst[k1]; k <= nlast[k1]; ++k) {
                ++kl;
                sum += cocc[kl + loopj] * cocc[kl + loopj];
            }
            aocc[kk] = sum;
        }
    }
    cutlim = 1.e-8;
    cutoff = std::max(cutlim, tiny * 10.0 * cutlim);
    flim = std::min(3.0, fref * 0.5);
    fref = 0.0;
    if (idiagg <= 5) {
        cutoff = cutlim;
    }
    sumt = 0.0;
    ijc = 0;
    tiny = 0.0;
    for (i = 1; i <= nvir; ++i) {
        loopi = ncvir[i];
        std::fill(latoms.begin(), latoms.end(), false);
        l = 0;
        for (j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j) {
            j1 = icvir[j];
            sum = 0.0;
            for (k = l + 1; k <= l + iorbs[j1]; ++k) sum += cvir[k + loopi] * cvir[k + loopi];
            l += iorbs[j1];
            avir[j1] = sum;
            latoms[icvir[j]] = true;
        }
        if (lijbo) {
            for (jj = nnce[i] + 1; jj <= nnce[i] + nce[i]; ++jj) {
                j1 = icvir[jj];
                for (jx = 1; jx <= iorbs[j1]; ++jx) ws[nfirst[j1] + jx - 1] = 0.0;
                kl = loopi;
                for (kk = nnce[i] + 1; kk <= nnce[i] + nce[i]; ++kk) {
                    k1 = icvir[kk];
                    kj = nijbo[k1][j1];
                    if (kj >= 0) {
                        if (avir[k1] * p[kj + 1] > cutoff) {
                            //  EXTRACT THE ATOM-ATOM INTERSECTION OF FAO
                            if (iorbs[k1] == 1 && iorbs[j1] == 1) {
                                ws[nfirst[j1]] += fao[kj + 1] * cvir[kl + 1];
                            } else if (k1 > j1) {
                                ii = kj;
                                for (i4 = 1; i4 <= iorbs[k1]; ++i4)
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        ++ii;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            } else if (k1 < j1) {
                                ii = kj;
                                for (jx = 1; jx <= iorbs[j1]; ++jx)
                                    for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                        ++ii;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            } else {
                                for (jx = 1; jx <= iorbs[j1]; ++jx)
                                    for (i4 = 1; i4 <= iorbs[j1]; ++i4) {
                                        if (i4 > jx) ii = kj + (i4 * (i4 - 1)) / 2 + jx;
                                        else ii = kj + (jx * (jx - 1)) / 2 + i4;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            }
                        }
                    }
                    kl += iorbs[k1];
                }
            }
        } else {
            for (jj = nnce[i] + 1; jj <= nnce[i] + nce[i]; ++jj) {
                j1 = icvir[jj];
                for (jx = 1; jx <= iorbs[j1]; ++jx) ws[nfirst[j1] + jx - 1] = 0.0;
                kl = loopi;
                for (kk = nnce[i] + 1; kk <= nnce[i] + nce[i]; ++kk) {
                    k1 = icvir[kk];
                    kj = ijbo(k1, j1);
                    if (kj >= 0) {
                        if (avir[k1] * p[kj + 1] > cutoff) {
                            if (iorbs[k1] == 1 && iorbs[j1] == 1) {
                                ws[nfirst[j1]] += fao[kj + 1] * cvir[kl + 1];
                            } else if (k1 > j1) {
                                ii = kj;
                                for (i4 = 1; i4 <= iorbs[k1]; ++i4)
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        ++ii;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            } else if (k1 < j1) {
                                ii = kj;
                                for (jx = 1; jx <= iorbs[j1]; ++jx)
                                    for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                        ++ii;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            } else {
                                for (jx = 1; jx <= iorbs[j1]; ++jx)
                                    for (i4 = 1; i4 <= iorbs[j1]; ++i4) {
                                        if (i4 > jx) ii = kj + (i4 * (i4 - 1)) / 2 + jx;
                                        else ii = kj + (jx * (jx - 1)) / 2 + i4;
                                        ws[nfirst[j1] + jx - 1] += fao[ii] * cvir[kl + i4];
                                    }
                            }
                        }
                    }
                    kl += iorbs[k1];
                }
            }
        }
        for (j = 1; j <= numat; ++j) {
            if (latoms[j]) {
                sum = 0.0;
                for (k = nfirst[j]; k <= nlast[j]; ++k) sum += ws[k] * ws[k];
                aov[j] = sum;
            } else {
                aov[j] = 0.0;
            }
        }
        //  EVALUATE THE VIRTUAL ENERGY LEVELS
        sum = 0.0;
        kl = loopi;
        for (kk = nnce[i] + 1; kk <= nnce[i] + nce[i]; ++kk) {
            k1 = icvir[kk];
            if (aov[k1] * avir[k1] > cutoff) {
                for (k = nfirst[k1]; k <= nlast[k1]; ++k) {
                    ++kl;
                    sum += ws[k] * cvir[kl];
                }
            } else {
                kl += iorbs[k1];
            }
        }
        eigv[nocc + i] = sum;
        //  EVALUATE THE OCCUPIED-VIRTUAL LMO INTERACTION ENERGIES
        if (idiagg <= 5 || (idiagg % 2) == 0) {
            nf = 0;
            i1 = icvir[nnce[i] + 1];
            i2 = (nce[i] > 1) ? icvir[nnce[i] + 2] : i1;
            if (lijbo) {
                for (j = 1; j <= nocc; ++j) {
                    if (ijc != nij) {
                        j1 = icocc[nncf[j] + 1];
                        j2 = (ncf[j] > 1) ? icocc[nncf[j] + 2] : j1;
                        i1j1 = nijbo[i1][j1];
                        i1j2 = nijbo[i1][j2];
                        i2j1 = nijbo[i2][j1];
                        i2j2 = nijbo[i2][j2];
                        if (i1j1 >= 0 || i1j2 >= 0 || i2j1 >= 0 || i2j2 >= 0) {
                            sum = 0.0;
                            if (i1j1 >= 0) sum += std::fabs(fao[i1j1 + 1]);
                            if (i1j2 >= 0) sum += std::fabs(fao[i1j2 + 1]);
                            if (i2j1 >= 0) sum += std::fabs(fao[i2j1 + 1]);
                            if (i2j2 >= 0) sum += std::fabs(fao[i2j2 + 1]);
                            if (sum >= flim) {
                                loopj = ncocc[j];
                                lij = false;
                                sum = 0.0;
                                kl = 0;
                                for (kk = nncf[j] + 1; kk <= nncf[j] + ncf[j]; ++kk) {
                                    k1 = icocc[kk];
                                    if (aocc[kk] * aov[k1] < cutoff || !latoms[k1]) {
                                        kl += iorbs[k1];
                                    } else {
                                        lij = true;
                                        for (k = nfirst[k1]; k <= nlast[k1]; ++k) {
                                            ++kl;
                                            sum += ws[k] * cocc[kl + loopj];
                                        }
                                    }
                                }
                                sumt += std::fabs(sum);
                                tiny = std::max(tiny, std::fabs(sum));
                                if (lij && std::fabs(sum) > oldlim) {
                                    ++nf;
                                    ++ijc;
                                    ifmo[1][ijc] = i;
                                    ifmo[2][ijc] = j;
                                    fmo[ijc] = sum;
                                }
                            }
                        }
                    }
                }
            } else {
                for (j = 1; j <= nocc; ++j) {
                    if (ijc != nij) {
                        j1 = icocc[nncf[j] + 1];
                        j2 = (ncf[j] > 1) ? icocc[nncf[j] + 2] : j1;
                        i1j1 = ijbo(i1, j1);
                        i1j2 = ijbo(i1, j2);
                        i2j1 = ijbo(i2, j1);
                        i2j2 = ijbo(i2, j2);
                        if (i1j1 >= 0 || i1j2 >= 0 || i2j1 >= 0 || i2j2 >= 0) {
                            sum = 0.0;
                            if (i1j1 >= 0) sum += std::fabs(fao[i1j1 + 1]);
                            if (i1j2 >= 0) sum += std::fabs(fao[i1j2 + 1]);
                            if (i2j1 >= 0) sum += std::fabs(fao[i2j1 + 1]);
                            if (i2j2 >= 0) sum += std::fabs(fao[i2j2 + 1]);
                            if (sum >= flim) {
                                loopj = ncocc[j];
                                lij = false;
                                sum = 0.0;
                                kl = 0;
                                for (kk = nncf[j] + 1; kk <= nncf[j] + ncf[j]; ++kk) {
                                    k1 = icocc[kk];
                                    if (aocc[kk] * aov[k1] < cutoff || !latoms[k1]) {
                                        kl += iorbs[k1];
                                    } else {
                                        lij = true;
                                        for (k = nfirst[k1]; k <= nlast[k1]; ++k) {
                                            ++kl;
                                            sum += ws[k] * cocc[kl + loopj];
                                        }
                                    }
                                }
                                sumt += std::fabs(sum);
                                tiny = std::max(tiny, std::fabs(sum));
                                if (lij && std::fabs(sum) > oldlim) {
                                    ++nf;
                                    ++ijc;
                                    ifmo[1][ijc] = i;
                                    ifmo[2][ijc] = j;
                                    fmo[ijc] = sum;
                                }
                            }
                        }
                    }
                }
            }
            nfmo[i] = nf;
        } else {
            ij0 = ijc + mydisp;
            for (jj = 1; jj <= nfmo[i]; ++jj) {
                if (ijc + mydisp == nij) break;
                j = ifmo[2][ij0 + jj];
                loopj = ncocc[j];
                sum = 0.0;
                kl = 0;
                for (kk = nncf[j] + 1; kk <= nncf[j] + ncf[j]; ++kk) {
                    k1 = icocc[kk];
                    if (aov[k1] * aocc[kk] < cutoff || !latoms[k1]) {
                        kl += iorbs[k1];
                    } else {
                        for (k = nfirst[k1]; k <= nlast[k1]; ++k) {
                            ++kl;
                            sum += ws[k] * cocc[kl + loopj];
                        }
                    }
                }
                sumt += std::fabs(sum);
                tiny = std::max(tiny, std::fabs(sum));
                ++ijc;
                fmo[ijc] = sum;
            }
        }
    }
    if (ijc == nij) {
        //   THERE WAS NOT ENOUGH STORAGE TO HOLD ALL THE INTEGRALS.
        //   THEREFORE, ON THE NEXT ITERATION, CALCULATE FEWER INTEGRALS.
        safety *= 2.0;
    } else {
        //  THERE IS ENOUGH STORAGE FOR ALL THE INTEGRALS.  IF NECESSARY,
        //  CALCULATE MORE INTEGRALS.
        safety = std::max(safety * 0.5, 1.0);
    }
    nij = ijc;
    if (idiagg > 2 && (idiagg % 4) != 0) {
        fref = tiny * tiny * tiny * tiny;
        if (times1) timer(" AFTER DIAGG1 IN ITER");
    } else {
        //  EVALUATE THE OCCUPIED ENERGY LEVELS
        for (i = 1; i <= nocc; ++i) {
            loopi = ncocc[i];
            l = 0;
            for (j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) {
                j1 = icocc[j];
                sum = 0.0;
                for (k = l + 1; k <= l + iorbs[j1]; ++k) sum += cocc[k + loopi] * cocc[k + loopi];
                l += iorbs[j1];
                avir[j1] = sum;
            }
            sum = 0.0;
            jl = loopi;
            if (lijbo) {
                for (j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) {
                    j1 = icocc[j];
                    kl = loopi;
                    for (k = nncf[i] + 1; k <= nncf[i] + ncf[i]; ++k) {
                        k1 = icocc[k];
                        k1j1 = nijbo[k1][j1];
                        if (k1j1 >= 0) {
                            if (avir[k1] * p[k1j1 + 1] * aocc[k] >= cutoff) {
                                if (k1 > j1) {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (i4 - 1) * iorbs[j1] + jx;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                } else if (k1 < j1) {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (jx - 1) * iorbs[k1] + i4;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                } else {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (j4 = 1; j4 <= jx; ++j4) {
                                            ii = k1j1 + (jx * (jx - 1)) / 2 + j4;
                                            sum1 += fao[ii] * cocc[kl + j4];
                                        }
                                        for (i4 = jx + 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (i4 * (i4 - 1)) / 2 + jx;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                }
                            }
                        }
                        kl += iorbs[k1];
                    }
                    jl += iorbs[j1];
                }
            } else {
                for (j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) {
                    j1 = icocc[j];
                    kl = loopi;
                    for (k = nncf[i] + 1; k <= nncf[i] + ncf[i]; ++k) {
                        k1 = icocc[k];
                        k1j1 = ijbo(k1, j1);
                        if (k1j1 >= 0) {
                            if (avir[k1] * p[k1j1 + 1] * aocc[k] >= cutoff) {
                                if (k1 > j1) {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (i4 - 1) * iorbs[j1] + jx;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                } else if (k1 < j1) {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (i4 = 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (jx - 1) * iorbs[k1] + i4;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                } else {
                                    for (jx = 1; jx <= iorbs[j1]; ++jx) {
                                        sum1 = 0.0;
                                        for (j4 = 1; j4 <= jx; ++j4) {
                                            ii = k1j1 + (jx * (jx - 1)) / 2 + j4;
                                            sum1 += fao[ii] * cocc[kl + j4];
                                        }
                                        for (i4 = jx + 1; i4 <= iorbs[k1]; ++i4) {
                                            ii = k1j1 + (i4 * (i4 - 1)) / 2 + jx;
                                            sum1 += fao[ii] * cocc[kl + i4];
                                        }
                                        sum += cocc[jl + jx] * sum1;
                                    }
                                }
                            }
                        }
                        kl += iorbs[k1];
                    }
                    jl += iorbs[j1];
                }
            }
            eigs[i] = sum;
        }
        oldlim = tiny * safety * 1.e-3;
        fref = tiny * tiny * tiny * tiny;
        if (times1) timer(" AFTER DIAGG1 IN ITER");
    }
    ovmax = tiny;
}
void diagg2(int nocc, int nvir, std::vector<double>& eigv,
            std::vector<int>& iused, std::vector<bool>& latoms, int nij,
            int idiagg, std::vector<double>& storei, std::vector<double>& storej) {
    // eigv is the eigs(nocc+1:) slice: element i (1..nvir) at eigv[nocc+i].
    bool retry;
    int i, ii, jur, l, ij, ilr, incv, iur, j, jlr, jncf, k, le, lf, lij,
        loopi, loopj, mie, mle, mlee, mlf, mlff, ncei, ncfj, nrej;
    double a, alpha, b, beta, biglim, c, d, e, sum;
    if (numcal != icalcn2) {
        icalcn2 = numcal;
        times2 = (keywrd.find(" TIMES") != std::string::npos);
        debug2 = (keywrd.find(" DIAGG2") != std::string::npos);
        //   IF THE SYSTEM IS A SOLID, THEN DAMP ROTATION OF VECTORS
        size_t pos_damp = keywrd.find(" DAMP");
        if (pos_damp != std::string::npos) {
            // Fortran: reada(keywrd, i+5) with 1-based i; 0-based pos = i-1,
            // so the 1-based read position is pos+6.
            const2 = reada(keywrd, (int)pos_damp + 6);
        } else if (id == 3) {
            const2 = 0.5;
        } else {
            a = 1.0;  //  If a transition metal, set DAMP to 0.5d0
            for (i = 1; i <= numat; ++i)
                if (!parameters_C::main_group[nat[i]]) a = 0.5;
            const2 = a;
        }
        //   EPS IS THE SMALLEST NUMBER WHICH, WHEN ADDED TO 1.D0, IS NOT EQUAL TO 1.D0
        epseta(eps2, eta2);
        //   INCREASE EPS TO ALLOW FOR A LOT OF ROUND-OFF
        bigeps = 50.0 * std::sqrt(eps2);
    }
    //  RETRY IS .TRUE. IF THE NUMBER OF REJECTED ANNIHILATIONS IS IN THE
    //  RANGE 1 TO 20 AND THE SAME ON TWO ITERATIONS.
    retry = (nrejct[0] == nrejct[1] && nrejct[0] != 0 && nrejct[0] < 20);
    if ((idiagg % 5) == 0 || idiagg <= 5) {
        tiny = -1.0;
        biglim = -1.0;
    } else {
        tiny = 0.01 * tiny;
        biglim = bigeps;
    }
    //   DO A CRUDE 2 BY 2 ROTATION TO "ELIMINATE" SIGNIFICANT ELEMENTS
    std::fill(iused.begin(), iused.end(), -1);
    std::fill(latoms.begin(), latoms.end(), false);
    if (debug2) {
        std::fprintf(stdout, "\n");
        std::fprintf(stdout, "            SIZE OF OCCUPIED ARRAYS IN DIAGG2\n\n");
        std::fprintf(stdout, "    LMO    NNCF     NCF   SPACE   NCOCC    SIZE   SPACE\n");
        for (i = 1; i <= nocc - 1; ++i) {
            l = ncocc[i];
            for (j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) l += iorbs[icocc[j]];
            std::fprintf(stdout, "%8d%8d%8d%8d%8d%8d%8d\n", i, nncf[i], ncf[i],
                         nncf[i + 1] - nncf[i] - ncf[i], ncocc[i], l - ncocc[i],
                         ncocc[i + 1] - l);
        }
        i = nocc;
        l = ncocc[i];
        for (j = nncf[i] + 1; j <= nncf[i] + ncf[i]; ++j) l += iorbs[icocc[j]];
        std::fprintf(stdout, "%8d%8d%8d%8d%8d%8d%8d\n", i, nncf[i], ncf[i],
                     icocc_dim - nncf[i] - ncf[i], ncocc[i], l - ncocc[i],
                     cocc_dim - l);
        std::fprintf(stdout, "\n");
        std::fprintf(stdout, "            SIZE OF VIRTUAL ARRAYS IN DIAGG2\n\n");
        std::fprintf(stdout, "    LMO    NNCE     NCE   SPACE   NCVIR    SIZE   SPACE\n");
        for (i = 1; i <= nvir - 1; ++i) {
            l = ncvir[i];
            for (j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j) l += iorbs[icvir[j]];
            std::fprintf(stdout, "%8d%8d%8d%8d%8d%8d%8d\n", i, nnce[i], nce[i],
                         nnce[i + 1] - nnce[i] - nce[i], ncvir[i], l - ncvir[i],
                         ncvir[i + 1] - l);
        }
        i = nvir;
        l = ncvir[i];
        for (j = nnce[i] + 1; j <= nnce[i] + nce[i]; ++j) l += iorbs[icvir[j]];
        std::fprintf(stdout, "%8d%8d%8d%8d%8d%8d%8d\n", i, nnce[i], nce[i],
                     icvir_dim - nnce[i] - nce[i], ncvir[i], l - ncvir[i],
                     cvir_dim - l);
    }
    sumb = 0.0;
    nrej = 0;
    lij = 0;
    for (ij = 1; ij <= nij; ++ij) {
        i = ifmo[1][ij];
        j = ifmo[2][ij];
        if (std::fabs(fmo[ij]) >= tiny) {
            c = fmo[ij] * const2;
            d = eigs[j] - eigv[nocc + i] - shift;
            if (std::fabs(c / d) >= biglim) {
                ncfj = ncf[j];
                ncei = nce[i];
                //  STORE LMOS FOR POSSIBLE REJECTION, IF LMOS EXPAND TOO MUCH.
                jlr = ncocc[j] + 1;
                if (j != nocc) {
                    jur = ncocc[j + 1];
                    jncf = nncf[j + 1];
                } else {
                    jur = cocc_dim;
                    jncf = icocc_dim;
                }
                jur = std::min(jlr + norbs - 1, jur);
                ilr = ncvir[i] + 1;
                if (i != nvir) {
                    iur = ncvir[i + 1];
                    incv = nnce[i + 1];
                } else {
                    iur = cvir_dim;
                    incv = icvir_dim;
                }
                iur = std::min(ilr + norbs - 1, iur);
                l = 0;
                for (k = jlr; k <= jur; ++k) { ++l; storej[l] = cocc[k]; }
                l = 0;
                for (k = ilr; k <= iur; ++k) { ++l; storei[l] = cvir[k]; }
                //   STORAGE DONE.
                ++lij;
                e = std::copysign(std::sqrt(4.0 * c * c + d * d), d);
                alpha = std::sqrt(0.5 * (1.0 + d / e));
                for (;;) {
                    beta = -std::copysign(std::sqrt(1.0 - alpha * alpha), c);
                    sumb += std::fabs(beta);
                    // IDENTIFY THE ATOMS IN THE OCCUPIED LMO.
                    mlf = 0;
                    for (lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf) {
                        ii = icocc[lf];
                        iused[ii] = mlf;
                        mlf += iorbs[ii];
                    }
                    loopi = ncvir[i];
                    loopj = ncocc[j];
                    mle = 0;
                    //      ROTATION OF PSEUDO-EIGENVECTORS
                    for (le = nnce[i] + 1; le <= nnce[i] + nce[i]; ++le) {
                        mie = icvir[le];
                        latoms[mie] = true;
                        mlff = iused[mie] + loopj;
                        if (iused[mie] >= 0) {
                            //  TWO BY TWO ROTATION OF ATOMS COMMON TO OCCUPIED LMO J AND VIRTUAL LMO I
                            for (mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee) {
                                ++mlff;
                                a = cocc[mlff];
                                b = cvir[mlee];
                                cocc[mlff] = alpha * a + beta * b;
                                cvir[mlee] = alpha * b - beta * a;
                            }
                        } else {
                            //   FILLED  LMO ATOM 'MIE' DOES NOT EXIST.
                            sum = 0.0;
                            for (mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee)
                                sum += (beta * cvir[mlee]) * (beta * cvir[mlee]);
                            if (sum > thresh) {
                                if (nncf[j] + ncf[j] >= jncf) goto label1000;
                                if (mlf + iorbs[mie] + loopj > jur) goto label1000;
                                //  YES, OCCUPIED LMO ATOM 'MIE' SHOULD EXIST
                                ++ncf[j];
                                icocc[nncf[j] + ncf[j]] = mie;
                                iused[mie] = mlf;
                                mlf += iorbs[mie];
                                //   PUT INTENSITY INTO OCCUPIED LMO ATOM 'MIE'
                                mlff = iused[mie] + loopj;
                                for (mlee = mle + 1 + loopi; mlee <= mle + iorbs[mie] + loopi; ++mlee) {
                                    ++mlff;
                                    cocc[mlff] = beta * cvir[mlee];
                                    cvir[mlee] = alpha * cvir[mlee];
                                }
                            }
                        }
                        mle += iorbs[mie];
                    }
                    //  NOW CHECK ALL ATOMS WHICH WERE IN THE OCCUPIED LMO
                    //  WHICH ARE NOT IN THE VIRTUAL LMO, TO SEE IF THEY
                    //  SHOULD BE IN THE VIRTUAL LMO.
                    for (lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf) {
                        ii = icocc[lf];
                        if (!latoms[ii]) {
                            sum = 0.0;
                            for (mlff = iused[ii] + loopj + 1; mlff <= iused[ii] + loopj + iorbs[ii]; ++mlff)
                                sum += (beta * cocc[mlff]) * (beta * cocc[mlff]);
                            if (sum > thresh) {
                                if (nnce[i] + nce[i] >= incv) goto label1000;
                                if (mle + iorbs[ii] + loopi > iur) goto label1000;
                                //  YES, VIRTUAL  LMO ATOM 'II' SHOULD EXIST
                                ++nce[i];
                                icvir[nnce[i] + nce[i]] = ii;
                                latoms[ii] = true;
                                //   PUT INTENSITY INTO VIRTUAL  LMO ATOM 'II'
                                mlff = iused[ii] + loopj;
                                for (mlee = mle + 1 + loopi; mlee <= mle + iorbs[ii] + loopi; ++mlee) {
                                    ++mlff;
                                    cvir[mlee] = -beta * cocc[mlff];
                                    cocc[mlff] = alpha * cocc[mlff];
                                }
                                mle += iorbs[ii];
                            }
                        }
                    }
                    break;
                label1000:
                    ++nrej;
                    //   THE ARRAY BOUNDS WERE GOING TO BE EXCEEDED.
                    //   TO PREVENT THIS, RESET THE LMOS.
                    l = 0;
                    for (k = jlr; k <= jur; ++k) { ++l; cocc[k] = storej[l]; }
                    l = 0;
                    for (k = ilr; k <= iur; ++k) { ++l; cvir[k] = storei[l]; }
                    ncf[j] = ncfj;
                    nce[i] = ncei;
                    for (k = 1; k <= numat; ++k) {
                        iused[k] = -1;
                        latoms[k] = false;
                    }
                    if (retry) {
                        //   HALF THE ROTATION ANGLE.
                        alpha = 0.5 * (alpha + 1.0);
                    } else {
                        goto next_ij;
                    }
                }
                //  RESET COUNTERS WHICH HAVE BEEN SET.
                for (le = nnce[i] + 1; le <= nnce[i] + nce[i]; ++le) {
                    mie = icvir[le];
                    latoms[mie] = false;
                }
                for (lf = nncf[j] + 1; lf <= nncf[j] + ncf[j]; ++lf)
                    iused[icocc[lf]] = -1;
            }
        }
    next_ij:;
    }
    nrejct[1] = nrejct[0];
    nrejct[0] = nrej;
    if (times2) timer(" AFTER DIAGG2 IN ITER");
}

void diagg(const std::vector<double>& fao, int nocc, int nvir, int idiagg,
           std::vector<double>& partp, int indi) {
    std::vector<double> storei(norbs), storej(norbs), ws(norbs), aov(numat + 1),
        avir(norbs);
    std::vector<double> aocc(std::max(1, icocc_dim));
    std::vector<int> iused(numat + 1);
    std::vector<bool> latoms(numat + 1);
    // allocate failures cannot occur with std::vector; memory errors throw.
    int nij = fmo_dim;
    //  In diagg1, all significant matrix elements are constructed.  These
    //  nij elements are stored in fmo, and their indices are stored in ifmo.
    diagg1(fao, nocc, nvir, eigs, ws, latoms, ifmo, fmo, fmo_dim, nij, idiagg,
           avir, aocc, aov);
    //  In diagg2, the significant matrix elements are annihilated by a two by
    //  two rotation.
    diagg2(nocc, nvir, eigs, iused, latoms, nij, idiagg, storei, storej);
    density_for_MOZYME(p, indi, nocc, partp);
}
