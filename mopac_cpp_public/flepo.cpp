// flepo.cpp — C++ translation of MOPAC 2016 "flepo.F90" (820 lines).
// BFGS (or DFP) variable-metric geometry optimiser. Minimises the heat of
// formation along search directions produced by the inverse-Hessian update;
// line searches by linmin; convergence by Herbert/gradient/xparam/function/
// Peters tests; periodic dumps via dfpsav. Fortran 1-based indexing kept.
#include "flepo.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

extern void compfg(const std::vector<double>& xparam, bool int_flag,
                   double& escf, bool fulscf, std::vector<double>& grad,
                   bool lgrad);
extern void linmin(double* xparam, double& alpha, double* pvect, int nvar,
                   double& funct, bool& okf, int& ic, double dott);
extern void supdot(double* s, const double* h, const double* g, int n,
                   int hoff = 0);
extern void dfpsav(double& totime, std::vector<double>& xparam,
                   std::vector<double>& gd, std::vector<double>& xlast,
                   double& funct1, std::vector<int>& mdfp,
                   std::vector<double>& xdfp);
extern void geout(int);
extern void mopend(const char* message);
extern void prttim(double tleft, double& tprt, char& txt);
extern double second(int n);
extern double reada(const std::string& s, int istart);
extern double ddot(int n, const double* x, int incx, const double* y,
                   int incy);
extern void dcopy(int n, const double* x, int incx, double* y, int incy);

void flepo_update_H(int nvar, const std::vector<double>& xvar,
                    const std::vector<double>& /*gvar*/,
                    const std::vector<double>& gg,
                    double sy, double yhy0, bool dfp,
                    std::vector<double>& hesinv) {
    double yhy = yhy0;
    if (!dfp) yhy = 1.0 + yhy / sy;
    int k = 0;
    for (int i = 1; i <= nvar; ++i) {
        double xvari = xvar[i] / sy;
        if (dfp) {
            double ggi = gg[i] / yhy0;
            for (int j = 1; j <= i; ++j)
                hesinv[k + j] += xvar[j] * xvari - gg[j] * ggi;
        } else {
            double ggi = gg[i] / sy;
            for (int j = 1; j <= i; ++j)
                hesinv[k + j] +=
                    -gg[j] * xvari - xvar[j] * ggi + yhy * xvar[j] * xvari;
        }
        k += i;
    }
}

void flepo(std::vector<double>& xparam, int nvar, double& funct1) {
    static int icalcn = 0;
    static double rst, sfact, dell, einc, del, const_s;
    static bool restrt, geook, saddle, minprt, dfp, thiel, print;
    static int igg1, nrst, ihdim, itry1, jcyc, lnstop, irepet, jnrst, ncount,
        maxcyc;
    static double rootv, delhof, tolerf, tolerg, tolrg, tolerx, drop, frepf,
        cncadd, absmin, alpha, pnorm, cycmx;

    emin = 0.0;
    double gnormr = 0.0, cos = 0.0, totime = 0.0, tlast = 0.0, tx1 = 0.0,
           tx2 = 0.0;
    bool limscf = true, resfil = false, lgrad = false, okf = false;
    int ic = 2;
    double beta = 0.0, smval = 0.0, dropn = 0.0;
    std::vector<double> gd, glast, xlast;  // (nvar+1), F90 flepo locals
    std::vector<double> xvar, gvar, xd, gg0, pvect0;  // scratch (nvar+1)
    if (icalcn != numcal) {
        hesinv.resize((nvar * (nvar + 1)) / 2 + 1, 0.0);
        rst = -1.0;
        nrst = 30;
        sfact = 1.5;
        dell = 0.01;
        einc = 0.3;
        igg1 = 3;
        del = dell;
        maxcyc = 100000;
        if (keywrd.find(" CYCLES") != std::string::npos)
            maxcyc = (int)std::lround(reada(
                keywrd, (int)keywrd.find(" CYCLES")));
        restrt = keywrd.find("RESTAR") != std::string::npos;
        thiel = keywrd.find("NOTHIE") == std::string::npos;
        geook = keywrd.find("GEO-OK") != std::string::npos;
        saddle = keywrd.find("SADDLE") != std::string::npos;
        minprt = !saddle;
        const_s = 1.0;
        dfp = keywrd.find("DFP") != std::string::npos;
        tolerg = 1.0;
        if (keywrd.find("PREC") != std::string::npos) tolerg = 0.2;
        if (keywrd.find("FORCE") != std::string::npos) tolerg = 0.1;
        if (keywrd.find("GNORM=") != std::string::npos) {
            rootv = 1.0;
            const_s = 1e-20;
            tolerg = reada(keywrd, (int)keywrd.find("GNORM="));
            if (keywrd.find(" LET") == std::string::npos && tolerg < 1e-2) {
                std::printf("  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01\n");
                tolerg = 1e-2;
            }
        } else {
            rootv = std::sqrt(nvar + 1e-5);
        }
        tolerx = 0.0001 * const_s;
        delhof = 0.0010 * const_s;
        tolerf = 0.002 * const_s;
        tolrg = tolerg;
        if (keywrd.find("PREC") != std::string::npos) {
            tolerx *= 0.01;
            delhof *= 0.01;
            tolerf *= 0.01;
            einc *= 0.01;
        }
        // therb computed per-call below
        tlast = tleft;
        tx2 = second(2);
        tleft = tleft - tx2 + time0;
        print = keywrd.find("FLEPO") != std::string::npos;
        drop = 1e15;
        frepf = 1e15;
        ihdim = (nvar * (nvar + 1)) / 2;
        cncadd = 1.0 / rootv;
        cncadd = std::min(0.15, cncadd);
        icalcn = numcal;
        if (restrt && nvar > 0) {
            std::vector<int> mdfp(10, 0);
            std::vector<double> xdfp(10, 0.0);
            gd.assign(nvar + 1, 0.0);
            xlast.assign(nvar + 1, 0.0);
            jnrst = 1;
            mdfp[9] = 0;
            mdfp[1] = jcyc; mdfp[2] = jnrst; mdfp[3] = ncount;
            mdfp[4] = lnstop;
            xdfp[1] = alpha; xdfp[2] = cos; xdfp[3] = pnorm; xdfp[4] = drop;
            xdfp[5] = del; xdfp[6] = frepf; xdfp[7] = cycmx; xdfp[8] = totime;
            dfpsav(totime, xparam, gd, xlast, funct1, mdfp, xdfp);
            if (moperr) return;
            jcyc = mdfp[1] - 1;
            jnrst = mdfp[2]; ncount = mdfp[3]; lnstop = mdfp[4];
            alpha = xdfp[1]; cos = xdfp[2]; pnorm = xdfp[3]; drop = xdfp[4];
            del = xdfp[5]; frepf = xdfp[6]; cycmx = xdfp[7]; totime = xdfp[8];
            int i = (int)(totime / 1000000.0);
            totime = totime - i * 1000000.0;
            time0 = time0 - totime;
            nscf = mdfp[5];
            std::printf(
                "\n\n          TOTAL TIME USED SO FAR:%13.2f SECONDS\n", totime);
            if (keywrd.find(" 1SCF") != std::string::npos) {
                last = 1;
                lgrad = keywrd.find(" GRAD") != std::string::npos;
                compfg(xparam, true, funct1, true, grad, lgrad);
                if (moperr) return;
                iflepo = 13;
                emin = 0.0;
                return;
            }
        }
    }

    // FIRST, INITIALIZE THE VARIABLES.
    int ireset = 0;
    double therb = 5.0 * tolerg * rootv;
    absmin = 1e6;
    itry1 = 0;
    jcyc = 0;
    lnstop = 1;
    irepet = 1;
    limscf = true;
    alpha = 1.0;
    pnorm = 1.0;
    jnrst = 0;
    cycmx = 0.0;
    cos = 0.0;
    totime = 0.0;
    ncount = 1;
    resfil = false;
    if (const_s > 1e-5 && saddle) {
        if (nvar > 0) gnorm = std::sqrt(ddot(nvar, grad.data(), 1,
                                             grad.data(), 1)) - 3.0;
        gnorm = std::min(10.0, gnorm);
        if (gnorm > 1.0) tolerg = tolrg * gnorm;
        std::printf(" GRADIENT CRITERION IN FLEPO =%10.3f\n", tolerg);
    }
    if (nvar == 1) {
        // Single variable: F90 goes to line 270 (linmin) directly, then the
        // nvar==1 block below; inlined here to avoid goto across scopes.
        pvect0.assign(2, 0.0);
        pvect0[1] = 0.01;
        alpha = 1.0;
        linmin(xparam.data(), alpha, pvect0.data(), nvar, funct1, okf, ic,
               dropn);
        if (moperr) return;
        std::printf(" ONLY ONE VARIABLE, THEREFORE ENERGY A MINIMUM\n");
        std::fflush(stdout);
        last = 1;
        lgrad = keywrd.find("GRAD") != std::string::npos;
        grad[1] = 0.0;
        compfg(xparam, true, funct1, true, grad, lgrad);
        if (moperr) return;
        iflepo = 14;
        emin = 0.0;
        tx2 = second(1);
        if (tx2 - time0 > tleft) {
            iflepo = -1;
            double totim = totime + second(1) - time0;
            std::vector<int> mdfp(10, 0);
            std::vector<double> xdfp(10, 0.0);
            mdfp[9] = 1;
            mdfp[5] = nscf;
            mdfp[1] = jcyc; mdfp[2] = jnrst; mdfp[3] = ncount; mdfp[4] = lnstop;
            xdfp[1] = alpha; xdfp[2] = cos; xdfp[3] = pnorm; xdfp[4] = drop;
            xdfp[5] = del; xdfp[6] = frepf; xdfp[7] = cycmx; xdfp[8] = totime;
            dfpsav(totim, xparam, gd, xlast, funct1, mdfp, xdfp);
            if (moperr) return;
        }
        return;
    }
    totime = 0.0;
    compfg(xparam, true, funct1, true, grad, true);
    iflepo = 16;
    if (nvar == 0) return;
    if (moperr) return;
    gd.assign(nvar + 1, 0.0);
    glast.assign(nvar + 1, 0.0);
    dcopy(nvar, grad.data(), 1, gd.data(), 1);
    if (nvar != 0) {
        gnorm = std::sqrt(ddot(nvar, grad.data(), 1, grad.data(), 1));
        gnormr = gnorm;
        if (lnstop != 1 && cos > rst &&
            (jnrst < nrst || !dfp) && restrt) {
            dcopy(nvar, gd.data(), 1, glast.data(), 1);
        } else {
            dcopy(nvar, grad.data(), 1, glast.data(), 1);
        }
    }
    if (gnorm < tolerg || nvar == 0) {
        iflepo = 2;
        compfg(xparam, true, funct1, true, grad, restrt);
        if (moperr) return;
        tx2 = second(1);
        emin = 0.0;
        return;
    }
    tx1 = second(2);
    tleft = tleft - tx1 + tx2;

    xvar.assign(nvar + 1, 0.0);
    gvar.assign(nvar + 1, 0.0);
    xd.assign(nvar + 1, 0.0);
    gg0.assign(nvar + 1, 0.0);
    pvect0.assign(nvar + 1, 0.0);
    gd.assign(nvar + 1, 0.0);
    glast.assign(nvar + 1, 0.0);
    xlast.assign(nvar + 1, 0.0);

    gnorm = std::sqrt(ddot(nvar, grad.data(), 1, grad.data(), 1));
    if (gnormr < 1e-10) gnormr = gnorm;
    int icyc = jcyc;
    goto L30;
L10:
    if (cos < rst) {
        for (int i = 1; i <= nvar; ++i) gd[i] = 0.5;
    }
L30:
    jcyc = jcyc + 1;
    jnrst = jnrst + 1;
    int i80 = 0;
    if (i80 == 1 || lnstop == 1 || cos <= rst || (jnrst >= nrst && dfp)) {
L50:
        for (int i = 1; i <= nvar; ++i) {
            double step = std::abs(grad[i]) * 0.0002;
            step = std::max(0.01, std::min(0.04, step));
            xd[i] = xparam[i] - (del >= 0 ? (grad[i] >= 0 ? del : -del)
                                          : (grad[i] >= 0 ? -del : del));
            // sign(del,grad(i)): sign of grad applied to magnitude del
            xd[i] = xparam[i] - std::copysign(del, grad[i]);
        }
        double funct2 = 0.0;
        compfg(xd, true, funct2, true, gd, true);
        if (moperr) return;
        if (!geook &&
            std::sqrt(ddot(nvar, gd.data(), 1, gd.data(), 1)) / gnorm > 10.0 &&
            gnorm > 20 && jcyc > 2) {
            del = del / 10.0;
            if (del >= 0.00005) {
                goto L50;
            } else {
                std::printf(" GRADIENTS OF OLD GEOMETRY, GNORM=%13.6f\n", gnorm);
                for (int i = 1; i <= nvar; ++i) std::printf("%12.6f", grad[i]);
                std::printf("\n");
                double gdnorm = std::sqrt(ddot(nvar, gd.data(), 1, gd.data(), 1));
                std::printf(" GRADIENTS OF NEW GEOMETRY, GNORM=%13.6f\n", gdnorm);
                for (int i = 1; i <= nvar; ++i) std::printf("%12.6f", gd[i]);
                std::printf(
                    "\n\n\n                    CALCULATION ABANDONED AT THIS POINT!\n");
                std::printf(
                    "\n\n          SMALL CHANGES IN INTERNAL\n"
                    "          COORDINATES ARE   \n"
                    "          CAUSING A LARGE CHANGE IN THE DISTANCE BETWEEN\n"
                    "          CHEMICALLY-BOUND ATOMS. THE GEOMETRY OPTIMIZATION\n"
                    "          PROCEDURE WOULD LIKELY PRODUCE INCORRECT RESULTS\n");
                geout(1);
                mopend("CALCULATION ABANDONED IN FLEPO");
                return;
            }
        }
        ncount = ncount + 1;
        for (int i = 1; i <= ihdim; ++i) hesinv[i] = 0.0;
        double sum = 0.0;
        int ii = 0, j = 0;
        if (funct2 < funct1) {
            for (int i = 1; i <= nvar; ++i) {
                ii = ii + i;
                double deltag = grad[i] - gd[i];
                double deltax = xparam[i] - xd[i];
                if (std::abs(deltag) < 0.001) deltag = 0.001;
                double ggd = std::max(1.0, std::abs(gd[i]));
                hesinv[ii] = std::min(0.1 / ggd, deltax / deltag);
                if (hesinv[ii] <= 0.0) continue;
                j = j + 1;
                sum = sum + hesinv[ii];
            }
        } else {
            for (int i = 1; i <= nvar; ++i) {
                ii = ii + i;
                double deltag = grad[i] - gd[i];
                double deltax = xparam[i] - xd[i];
                if (std::abs(deltag) < 0.001) deltag = 0.001;
                double ggd = std::max(1.0, std::abs(grad[i]));
                hesinv[ii] = std::min(0.1 / ggd, deltax / deltag);
                if (hesinv[ii] <= 0.0) continue;
                j = j + 1;
                sum = sum + hesinv[ii];
            }
        }
        if (j != 0) {
            sum = sum / j;
            ii = 0;
            for (int i = 1; i <= nvar; ++i) {
                ii = ii + i;
                if (hesinv[ii] >= 0) continue;
                hesinv[ii] = sum;
            }
        }
        jnrst = 0;
        if (jcyc < 2) cosine = 1.0;
        if (funct2 >= funct1) {
            if (print)
                std::printf(" FUNCTION VALUE=%13.7f  WILL NOT BE REPLACED BY VALUE=%13.7f\n"
                            "          CALCULATED BY RESTART PROCEDURE\n\n",
                            funct1, funct2);
            cosine = 1.0;
        } else {
            if (print)
                std::printf(" FUNCTION VALUE=%13.7f IS BEING REPLACED BY VALUE=%13.7f\n"
                            "          FOUND IN RESTART PROCEDURE\n"
                            "          THE CORRESPONDING X VALUES AND GRADIENTS ARE ALSO BEING REPLACED\n\n",
                            funct1, funct2);
            funct1 = funct2;
            for (int i = 1; i <= nvar; ++i) {
                xparam[i] = xd[i];
                grad[i] = gd[i];
            }
            gnorm = std::sqrt(ddot(nvar, grad.data(), 1, grad.data(), 1));
            if (gnormr < 1e-10) gnormr = gnorm;
        }
    } else {
        // UPDATE VARIABLE-METRIC MATRIX
        for (int i = 1; i <= nvar; ++i) {
            xvar[i] = xparam[i] - xlast[i];
            gvar[i] = grad[i] - glast[i];
        }
        supdot(gg0.data(), hesinv.data(), gvar.data(), nvar);
        double yhy = ddot(nvar, gg0.data(), 1, gvar.data(), 1);
        double sy = ddot(nvar, xvar.data(), 1, gvar.data(), 1);
        int k = 0;
        if (dfp) {
            for (int i = 1; i <= nvar; ++i) {
                double xvari = xvar[i] / sy;
                double ggi = gg0[i] / yhy;
                for (int j = 1; j <= i; ++j)
                    hesinv[k + j] += xvar[j] * xvari - gg0[j] * ggi;
                k = i + k;
            }
        } else {
            yhy = 1.0 + yhy / sy;
            for (int i = 1; i <= nvar; ++i) {
                double xvari = xvar[i] / sy;
                double ggi = gg0[i] / sy;
                for (int j = 1; j <= i; ++j)
                    hesinv[k + j] +=
                        -gg0[j] * xvari - xvar[j] * ggi + yhy * xvar[j] * xvari;
                k = i + k;
            }
        }
    }
    // ESTABLISH NEW SEARCH DIRECTION
    double pnlast = pnorm;
    supdot(pvect0.data(), hesinv.data(), grad.data(), nvar);
    pnorm = std::sqrt(ddot(nvar, pvect0.data(), 1, pvect0.data(), 1));
    if (pnorm > 1.5 * pnlast) {
        for (int i = 1; i <= nvar; ++i) pvect0[i] = pvect0[i] * 1.5 * pnlast / pnorm;
        pnorm = 1.5 * pnlast;
    }
    double dott = -ddot(nvar, pvect0.data(), 1, grad.data(), 1);
    for (int i = 1; i <= nvar; ++i) pvect0[i] = -pvect0[i];
    cos = -dott / (pnorm * gnorm);
    if (print) {
        std::printf(" AT THE BEGINNING OF CYCLE%5d  THE FUNCTION VALUE IS %13.6f\n"
                    "  THE CURRENT POINT IS ...\n", jcyc, funct1);
        std::printf("  GRADIENT NORM = %10.4f\n  ANGLE COSINE =%10.4f\n",
                    gnorm, cos);
        int nto6 = (nvar - 1) / 6 + 1;
        int iinc1 = -5;
        for (int i = 1; i <= nto6; ++i) {
            std::printf("\n");
            iinc1 = iinc1 + 6;
            int iinc2 = std::min(iinc1 + 5, nvar);
            std::printf("     I");
            for (int j = iinc1; j <= iinc2; ++j) std::printf("%9d", j);
            std::printf("\n");
            std::printf("  XPARAM(I)");
            for (int j = iinc1; j <= iinc2; ++j) std::printf("%9.4f", xparam[j]);
            std::printf("\n");
            std::printf("  GRAD  (I)");
            for (int j = iinc1; j <= iinc2; ++j) std::printf("%10.4f", grad[j]);
            std::printf("\n");
            std::printf("  PVECT (I)");
            for (int j = iinc1; j <= iinc2; ++j) std::printf("%10.6f", pvect0[j]);
            std::printf("\n");
        }
    }
    lnstop = 0;
    alpha = alpha * pnlast / pnorm;
    for (int i = 1; i <= nvar; ++i) {
        glast[i] = grad[i];
        xlast[i] = xparam[i];
    }
    if (jnrst == 0) alpha = 1.0;
    drop = std::abs(alpha * dott);
    if (print) std::printf(" -ALPHA.P.G =%18.6f\n", drop);
    if (gnorm <= therb && jnrst != 0 && drop < delhof) {
        if (minprt)
            std::printf("\n\n          HERBERTS TEST SATISFIED - GEOMETRY OPTIMIZED\n");
        last = 1;
        compfg(xparam, true, funct1, true, grad, false);
        if (moperr) return;
        iflepo = 3;
        time0 = time0 - totime;
        emin = 0.0;
        tx2 = second(1);
        return;
    }
    beta = alpha;
    smval = funct1;
    dropn = -std::abs(drop / alpha);
    okf = false;
    ic = 2;
L270:
    linmin(xparam.data(), alpha, pvect0.data(), nvar, funct1, okf, ic, dropn);
    if (moperr) return;
    if (nvar == 1) {
        std::printf(" ONLY ONE VARIABLE, THEREFORE ENERGY A MINIMUM\n");
        last = 1;
        lgrad = keywrd.find("GRAD") != std::string::npos;
        grad[1] = 0.0;
        compfg(xparam, true, funct1, true, grad, lgrad);
        if (moperr) return;
        iflepo = 14;
        emin = 0.0;
        tx2 = second(1);
        if (tx2 - time0 > tleft) {
            iflepo = -1;
            double totim = totime + second(1) - time0;
            std::vector<int> mdfp(10, 0);
            std::vector<double> xdfp(10, 0.0);
            mdfp[9] = 1;
            mdfp[5] = nscf;
            mdfp[1] = jcyc; mdfp[2] = jnrst; mdfp[3] = ncount; mdfp[4] = lnstop;
            xdfp[1] = alpha; xdfp[2] = cos; xdfp[3] = pnorm; xdfp[4] = drop;
            xdfp[5] = del; xdfp[6] = frepf; xdfp[7] = cycmx; xdfp[8] = totime;
            dfpsav(totim, xparam, gd, xlast, funct1, mdfp, xdfp);
            if (moperr) return;
        }
        return;
    }
    if (ireset > 10 || (gnorm < 40.0 && gnorm / gnormr < 0.33)) {
        ireset = 0;
        gnormr = 0.0;
        for (int i = 1; i <= nvar; ++i) grad[i] = 0.0;
    }
    ireset = ireset + 1;
    if (thiel) {
        double sum = 0.0;
        compfg(xparam, ic != 1, sum, true, grad, true);
        if (moperr) return;
    } else {
        compfg(xparam, true, funct1, true, grad, true);
        if (moperr) return;
    }
    gnorm = std::sqrt(ddot(nvar, grad.data(), 1, grad.data(), 1));
    if (gnormr < 1e-10) gnormr = gnorm;
    ncount = ncount + 1;
    if (!okf) {
        lnstop = 1;
        if (minprt)
            std::printf(
                "\n                    NO POINT LOWER IN ENERGY THAN THE STARTING POINT\n"
                "                    COULD BE FOUND IN THE LINE MINIMIZATION\n");
        funct1 = smval;
        alpha = beta;
        for (int i = 1; i <= nvar; ++i) {
            grad[i] = glast[i];
            xparam[i] = xlast[i];
        }
        if (jnrst == 0) {
            std::printf(
                "\n\n                    SINCE COS WAS JUST RESET, THE SEARCH IS BEING ENDED\n");
            last = 1;
            compfg(xparam, true, funct1, true, grad, true);
            if (moperr) return;
            iflepo = 4;
            time0 = time0 - totime;
            tx2 = second(1);
            emin = 0.0;
            return;
        }
        if (print)
            std::printf("                     COS WILL BE RESET AND ANOTHER ATTEMPT MADE\n");
        cos = 0.0;
        // Hoisted declarations (C2362: goto L430 must not skip initialization).
        {
          double xn = 0.0, tx = 0.0, tf = 0.0;
          (void)xn; (void)tx; (void)tf;
        }
        goto L430;
    }
    double xn = std::sqrt(ddot(nvar, xparam.data(), 1, xparam.data(), 1));
    double tx = std::abs(alpha * pnorm);
    if (xn != 0.0) tx = tx / xn;
    double tf = std::abs(smval - funct1);
    if (absmin - smval < 1e-7) {
        itry1 = itry1 + 1;
        if (itry1 > 10) {
            std::printf("\n\n HEAT OF FORMATION IS ESSENTIALLY STATIONARY\n");
            goto L420;
        }
    } else {
        itry1 = 0;
        absmin = smval;
    }
    if (print)
        std::printf(
            "\n           NUMBER OF COUNTS =%6d         COS    =%11.4f\n"
            "  ABSOLUTE  CHANGE IN X     =%13.6f  ALPHA  =%11.4f\n"
            "  PREDICTED CHANGE IN F     =  %11.4f  ACTUAL =  %11.4f\n"
            "  GRADIENT NORM             =  %11.4f\n\n",
            ncount, cos, tx * xn, alpha, -drop, -tf, gnorm);
    if (tx <= tolerx) {
        if (minprt) std::printf(" TEST ON X SATISFIED\n");
        goto L350;
    }
    if (tf <= tolerf) {
        if (minprt) std::printf(" HEAT OF FORMATION TEST SATISFIED\n");
        goto L350;
    }
    if (gnorm <= tolerg * rootv) {
        if (minprt) std::printf(" TEST ON GRADIENT SATISFIED\n");
        goto L350;
    }
    goto L430;
L350:
    for (int i = 1; i <= nvar; ++i) {
        if (std::abs(grad[i]) <= tolerg) continue;
        irepet = irepet + 1;
        if (irepet <= 1) {
            frepf = funct1;
            cos = 0.0;
        }
        if (minprt)
            std::printf("                    HOWEVER, A COMPONENT OF GRADIENT IS LARGER THAN%6.2f\n\n",
                        tolerg);
        if (std::abs(funct1 - frepf) > einc) irepet = 0;
        if (irepet > igg1) {
            std::printf(
                "          THERE HAVE BEEN%2d ATTEMPTS TO REDUCE THE GRADIENT.\n"
                "          DURING THESE ATTEMPTS THE ENERGY DROPPED BY LESS THAN%4.1f KCAL/MOLE\n"
                "          FURTHER CALCULATION IS NOT JUSTIFIED AT THIS TIME.\n",
                igg1, einc);
            if (keywrd.find("PREC") == std::string::npos)
                std::printf("          TO CONTINUE, START AGAIN WITH THE WORD \"PRECISE\"\n");
            last = 1;
            compfg(xparam, true, funct1, true, grad, false);
            if (moperr) return;
            iflepo = 8;
            time0 = time0 - totime;
            tx2 = second(1);
            emin = 0.0;
            return;
        } else {
            goto L430;
        }
    }
    if (minprt) std::printf("PETERS TEST SATISFIED\n");
L420:
    last = 1;
    compfg(xparam, true, funct1, true, grad, false);
    if (moperr) return;
    iflepo = 6;
    time0 = time0 - totime;
    tx2 = second(1);
    emin = 0.0;
    return;
L430:
    double bsmvf = std::abs(smval - funct1);
    if (bsmvf > 10.0) cos = 0.0;
    del = 0.002;
    if (bsmvf > 1.0) del = dell / 2.0;
    if (bsmvf > 5.0) del = dell;
    tx2 = second(2);
    double tcycle = tx2 - tx1;
    tx1 = tx2;
    if (tcycle < 100000.0) cycmx = std::max(cycmx, tcycle);
    tleft = tleft - tcycle;
    if (tleft < 0) tleft = -0.1;
    if (tcycle > 1e5) tcycle = 0.0;
    if (tlast - tleft > tdump) {
        double totim = totime + second(1) - time0;
        tlast = tleft;
        std::vector<int> mdfp(10, 0);
        std::vector<double> xdfp(10, 0.0);
        mdfp[9] = 2;
        resfil = true;
        mdfp[5] = nscf;
        mdfp[1] = jcyc; mdfp[2] = jnrst; mdfp[3] = ncount; mdfp[4] = lnstop;
        xdfp[1] = alpha; xdfp[2] = cos; xdfp[3] = pnorm; xdfp[4] = drop;
        xdfp[5] = del; xdfp[6] = frepf; xdfp[7] = cycmx; xdfp[8] = totime;
        dfpsav(totim, xparam, gd, xlast, funct1, mdfp, xdfp);
        if (moperr) return;
    }
    double tprt = 0.0;
    char txt = 'S';
    prttim(tleft, tprt, txt);
    if (resfil) {
        if (minprt)
            std::printf(" RESTART FILE WRITTEN,      TIME LEFT:%6.2f%c  GRAD.:%10.3f HEAT:%14.7f\n",
                        tprt, txt, std::min(gnorm, 999999.999), funct1);
        if (chanel_C::log) std::fprintf(stderr, " CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%c  GRAD.:%10.3f HEAT:%14.7f\n",
                              jcyc, std::min(tcycle, 9999.99), tprt, txt,
                              std::min(gnorm, 999999.999), funct1);
        resfil = false;
    } else {
        if (minprt)
            std::printf(" CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%c  GRAD.:%10.3f HEAT:%14.7f\n",
                        jcyc, std::min(tcycle, 9999.99), tprt, txt,
                        std::min(gnorm, 999999.999), funct1);
        if (chanel_C::log) std::fprintf(stderr, " CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%c  GRAD.:%10.3f HEAT:%14.7f\n",
                              jcyc, std::min(tcycle, 9999.99), tprt, txt,
                              std::min(gnorm, 999999.999), funct1);
    }
    if (tleft > sfact * cycmx && jcyc - icyc < maxcyc) goto L10;
    std::printf("                    THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE\n"
                "                    NOW GOING TO FINAL\n");
    double totim = totime + second(1) - time0;
    std::vector<int> mdfp(10, 0);
    std::vector<double> xdfp(10, 0.0);
    mdfp[9] = 1;
    mdfp[5] = nscf;
    mdfp[1] = jcyc; mdfp[2] = jnrst; mdfp[3] = ncount; mdfp[4] = lnstop;
    xdfp[1] = alpha; xdfp[2] = cos; xdfp[3] = pnorm; xdfp[4] = drop;
    xdfp[5] = del; xdfp[6] = frepf; xdfp[7] = cycmx; xdfp[8] = totime;
    dfpsav(totim, xparam, gd, xlast, funct1, mdfp, xdfp);
    if (moperr) return;
    iflepo = -1;
    return;
}
