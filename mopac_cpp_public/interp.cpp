// interp.cpp
#include "interp.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

#include "chanel_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "rsp.h"
#include "schmib.h"
#include "schmit.h"

using funcon_C::fpc_9;
using funcon_C::pi;
using molkst_C::keywrd;
using molkst_C::numcal;
using molkst_C::norbs;

void spline(double* x, double* f, double* df, double xhigh, double xlow,
            double& xmin, int npnts) {
    const double close=1e-8, big=500.0, huge=1e10, ustep=1.0, dstep=2.0;
    double fmin=f[npnts], dfmin=df[npnts];
    int n1=npnts-1;
    int k = n1;
    // move nth point to sorted position
    double xm=x[npnts], fm=f[npnts], dfm=df[npnts];
    while (k>=1 && x[k] > xm) { x[k+1]=x[k]; f[k+1]=f[k]; df[k+1]=df[k]; --k; }
    x[k+1]=xm; f[k+1]=fm; df[k+1]=dfm;
    if (df[1] > 0) {}
    double step = (df[1] > 0) ? dstep : ustep;
    double xstart = x[1] - step*(x[2]-x[1]);
    xstart = std::max(xstart, xlow);
    step = (df[npnts] > 0) ? ustep : dstep;
    double xstop = x[npnts] + step*(x[npnts]-x[n1]);
    xstop = std::min(xstop, xhigh);
    xmin = xm; fmin = fm; dfmin = dfm;
    for (int kk=1; kk<=n1; ++kk) {
        if (f[kk] < fmin) { xmin=x[kk]; fmin=f[kk]; dfmin=df[kk]; }
        double dx = x[kk+1]-x[kk];
        if (dx <= close) continue;
        double x1=0.0, x2=dx;
        if (kk==1) x1=xstart-x[1];
        if (kk==n1) x2=xstop-x[n1];
        double dum=(f[kk+1]-f[kk])/dx;
        double a=(df[kk]+df[kk+1]-2*dum)/(dx*dx);
        double b=(3*dum-2*df[kk]-df[kk+1])/dx;
        double c=df[kk];
        double bb=b*b, ac3=3*a*c;
        double xk;
        if (bb < ac3) continue;
        if (b<=0) { if (std::fabs(b) > huge*std::fabs(a)) continue; }
        else { if (bb > big*std::fabs(ac3)) { double r=ac3/bb; xk=-(((0.039063*r+0.0625)*r+0.125)*r+0.5)*c/b; goto check; } }
        xk = (-b + std::sqrt(bb-ac3))/(3*a);
check:
        if (xk<x1 || xk>x2) continue;
        double fm2 = ((a*xk+b)*xk+c)*xk+f[kk];
        if (fm2 <= fmin) { xmin=xk+x[kk]; fmin=fm2; dfmin=(3*a*xk+2*b)*xk+c; }
    }
}

void interp(int np, int nq, int& mode, double e, double* fp, double* cp,
            double* theta, double* vec_interp, double* fock_interp,
            double* p_interp, double* h_interp, double* vecl,
            double& eold_ref) {
    static int icalcn = 0;
    static bool debug = false;
    static double emin = 0.0;
    static double xold = 0.0;
    static int npnts = 0;
    const double zero = 0.0, ff = 0.9;

    if (icalcn != numcal) {
        debug = (keywrd.find("INTERP") != std::string::npos);
        icalcn = numcal;
        emin = 0.0;
        npnts = 0;
    }

    const int minpq = std::min(np, nq);
    const int np1 = np + 1;
    const int np2 = std::max(1, np / 2);
    int ipoint = 0, ii = 0;
    double xnow = 0.0, xhigh = 0.0, xlow = 0.0;

    if (mode != 2) {
        // (MODE=1 OR 3 ENTRY) Transform Fock matrix to current MO basis.
        for (int i = 1; i <= norbs; ++i) {
            const int i1 = i + 1;
            for (int j = 1; j <= nq; ++j) {
                double dum = zero;
                for (int r = 1; r <= i; ++r)
                    dum += fp[ii + r - 1] * cp[(j + np - 1) * norbs + (r - 1)];
                if (i != norbs) {
                    int ik = ii + i + i;
                    for (int k = i1; k <= norbs; ++k) {
                        dum += fp[ik - 1] * cp[(j + np - 1) * norbs + (k - 1)];
                        ik += k;
                    }
                }
                p_interp[(j - 1) * norbs + (i - 1)] = dum;
            }
            ii += i;
        }
        for (int i = 1; i <= np; ++i) {
            for (int j = 1; j <= nq; ++j) {
                double dum = zero;
                for (int r = 1; r <= norbs; ++r)
                    dum += cp[(i - 1) * norbs + (r - 1)] *
                           p_interp[(j - 1) * norbs + (r - 1)];
                fock_interp[(j - 1) * norbs + (i - 1)] = dum;
            }
        }
        if (mode != 3) {
            // Current point becomes old point (MODE=1 entry).
            for (int c = 0; c < norbs * norbs; ++c) vec_interp[c] = cp[c];
            eold_ref = e;
            xold = 1.0;
            mode = 2;
            return;
        }
        // (MODE=3 ENTRY) F corresponds to current point.
        ++npnts;
        if (debug)
            std::printf("   INTERPOLATED ENERGY:%13.6f\n", e * fpc_9);
        ipoint = npnts;
    } else {
        // (MODE=2 ENTRY) Calculate theta, and U, V, W matrices.
        int j1 = 1;
        for (int i = 1; i <= norbs; ++i) {
            if (i == np1) j1 = np1;
            for (int j = j1; j <= norbs; ++j) {
                double s = zero;
                for (int r = 1; r <= norbs; ++r)
                    s += cp[(i - 1) * norbs + (r - 1)] *
                         vec_interp[(j - 1) * norbs + (r - 1)];
                p_interp[(j - 1) * norbs + (i - 1)] = s;
            }
        }
        // U = CP(dagger)*VEC is now in P array.
        int ij = 0;
        for (int i = 1; i <= np; ++i) {
            for (int j = 1; j <= i; ++j) {
                ++ij;
                double s = zero;
                for (int r = np1; r <= norbs; ++r)
                    s += p_interp[(r - 1) * norbs + (i - 1)] *
                         p_interp[(r - 1) * norbs + (j - 1)];
                h_interp[ij - 1] = s;
            }
        }
        rsp(h_interp, np, theta, vecl);
        // vec_interp(np:1:-1, i) = vecl(np+il : il+1 : -1), il=(i-1)*np
        for (int i = 1; i <= np; ++i) {
            const int il = i * np - np;
            for (int k = 1; k <= np; ++k)
                vec_interp[(i - 1) * norbs + (k - 1)] =
                    vecl[(i - 1) * np + (np - k)];
        }
        for (int i = 1; i <= np2; ++i) {
            double dum = theta[np1 - i - 1];
            theta[np1 - i - 1] = theta[i - 1];
            theta[i - 1] = dum;
            for (int j = 1; j <= np; ++j) {
                dum = vec_interp[(np1 - i - 1) * norbs + (j - 1)];
                vec_interp[(np1 - i - 1) * norbs + (j - 1)] =
                    vec_interp[(i - 1) * norbs + (j - 1)];
                vec_interp[(i - 1) * norbs + (j - 1)] = dum;
            }
        }
        for (int i = 1; i <= minpq; ++i) {
            theta[i - 1] = std::max(theta[i - 1], zero);
            theta[i - 1] = std::min(theta[i - 1], 1.0);
            theta[i - 1] = std::asin(std::sqrt(theta[i - 1]));
        }
        // Now compute WQ.
        for (int i = 1; i <= nq; ++i) {
            for (int j = 1; j <= minpq; ++j) {
                double s = zero;
                for (int r = 1; r <= np; ++r)
                    s += p_interp[(np + i - 1) * norbs + (r - 1)] *
                         vec_interp[(j - 1) * norbs + (r - 1)];
                vec_interp[(np + j - 1) * norbs + (i - 1)] = s;
            }
        }
        schmit(&vec_interp[(np1 - 1) * norbs], nq, norbs);
        // Transpose NP by NP block of U stored in P.
        for (int i = 1; i <= np; ++i) {
            for (int j = 1; j <= i; ++j) {
                double dum = p_interp[(j - 1) * norbs + (i - 1)];
                p_interp[(j - 1) * norbs + (i - 1)] =
                    p_interp[(i - 1) * norbs + (j - 1)];
                p_interp[(i - 1) * norbs + (j - 1)] = dum;
            }
        }
        // Calculate WP matrix, hold in first NP columns of P.
        for (int i = 1; i <= np; ++i) {
            for (int r = 1; r <= np; ++r)
                h_interp[r - 1] = p_interp[(r - 1) * norbs + (i - 1)];
            for (int j = 1; j <= np; ++j) {
                double s = zero;
                for (int r = 1; r <= np; ++r)
                    s += h_interp[r - 1] * vec_interp[(j - 1) * norbs + (r - 1)];
                p_interp[(j - 1) * norbs + (i - 1)] = s;
            }
        }
        schmib(p_interp, np, norbs);
        // Calculate VQ matrix, hold in last NQ columns of P matrix.
        for (int i = 1; i <= nq; ++i) {
            for (int r = 1; r <= nq; ++r)
                h_interp[r - 1] =
                    p_interp[(np + r - 1) * norbs + (np + i - 1)];
            for (int j = np1; j <= norbs; ++j) {
                double s = zero;
                for (int r = 1; r <= nq; ++r)
                    s += h_interp[r - 1] * vec_interp[(j - 1) * norbs + (r - 1)];
                p_interp[(j - 1) * norbs + (i - 1)] = s;
            }
        }
        schmib(&p_interp[(np1 - 1) * norbs], nq, norbs);
        // Calculate (DE/DX) at old point.
        double dedx = zero;
        for (int i = 1; i <= np; ++i) {
            for (int j = 1; j <= nq; ++j) {
                double dum = zero;
                for (int r = 1; r <= minpq; ++r)
                    dum += theta[r - 1] * p_interp[(r - 1) * norbs + (i - 1)] *
                           vec_interp[(np + r - 1) * norbs + (j - 1)];
                dedx += dum * fock_interp[(j - 1) * norbs + (i - 1)];
            }
        }
        const double deold = -4.0 * dedx;
        double x[13], f[13], df[13];
        x[2] = xold;
        f[2] = eold_ref;
        df[2] = deold;
        // Move VP out of vec_interp into first NP columns of p_interp.
        for (int c = 0; c < np * norbs; ++c) p_interp[c] = vec_interp[c];
        int k1 = 0, k2 = np;
        for (int j = 1; j <= norbs; ++j) {
            if (j == np1) { k1 = np; k2 = nq; }
            for (int i = 1; i <= norbs; ++i) {
                double dum = zero;
                for (int r = 1; r <= k2; ++r)
                    dum += cp[(k1 + r - 1) * norbs + (i - 1)] *
                           p_interp[(j - 1) * norbs + (r - 1)];
                vec_interp[(j - 1) * norbs + (i - 1)] = dum;
            }
        }
        // vec_interp(dagger)*FP*vec_interp; store off-diagonal block.
        ii = 0;
        for (int i = 1; i <= norbs; ++i) {
            const int i1 = i + 1;
            for (int j = 1; j <= nq; ++j) {
                double dum = zero;
                for (int r = 1; r <= i; ++r)
                    dum += fp[ii + r - 1] * vec_interp[(j + np - 1) * norbs + (r - 1)];
                if (i != norbs) {
                    int ik = ii + i + i;
                    for (int k = i1; k <= norbs; ++k) {
                        dum += fp[ik - 1] * vec_interp[(j + np - 1) * norbs + (k - 1)];
                        ik += k;
                    }
                }
                p_interp[(j - 1) * norbs + (i - 1)] = dum;
            }
            ii += i;
        }
        for (int i = 1; i <= np; ++i) {
            for (int j = 1; j <= nq; ++j) {
                double dum = zero;
                for (int r = 1; r <= norbs; ++r)
                    dum += vec_interp[(i - 1) * norbs + (r - 1)] *
                           p_interp[(j - 1) * norbs + (r - 1)];
                fock_interp[(j - 1) * norbs + (i - 1)] = dum;
            }
        }
        // Set limits on range of 1-D search.
        npnts = 2;
        ipoint = 1;
        xnow = zero;
        xhigh = pi / (2.0 * theta[0]);
        xlow = -0.5 * xhigh;
        // Current-point derivative and spline fit (mode=2 branch tail).
        dedx = zero;
        for (int k = 1; k <= minpq; ++k)
            dedx += theta[k - 1] * fock_interp[(k - 1) * norbs + (k - 1)];
        {
            const double denow2 = -4.0 * dedx;
            x[1] = xnow;
            f[1] = e;
            df[1] = denow2;
        }
        double xmin = 0.0;
        spline(x, f, df, xhigh, xlow, xmin, npnts);
        if (eold_ref - e <= ff * (eold_ref - emin) && ipoint <= 10) {
            const double xnow2 = xmin;
            for (int c = 0; c < norbs * norbs; ++c) cp[c] = vec_interp[c];
            for (int k = 1; k <= minpq; ++k) {
                const double ck = std::cos(xnow2 * theta[k - 1]);
                const double sk = std::sin(xnow2 * theta[k - 1]);
                if (debug)
                    std::printf(" ROTATION ANGLE:%12.4f\n", sk * 57.29578);
                for (int r = 1; r <= norbs; ++r) {
                    cp[(k - 1) * norbs + (r - 1)] =
                        ck * vec_interp[(k - 1) * norbs + (r - 1)] -
                        sk * vec_interp[(np + k - 1) * norbs + (r - 1)];
                    cp[(np + k - 1) * norbs + (r - 1)] =
                        sk * vec_interp[(k - 1) * norbs + (r - 1)] +
                        ck * vec_interp[(np + k - 1) * norbs + (r - 1)];
                }
            }
            mode = 3;
            return;
        }
        if (mode != 2) {
            for (int c = 0; c < norbs * norbs; ++c) vec_interp[c] = cp[c];
            mode = 2;
        }
        eold_ref = e;
        if (npnts <= 200) return;
        std::printf("  K     X(K)        F(K)         DF(K)\n");
        for (int k = 1; k <= npnts; ++k)
            std::printf("%3d%10.5f%15.10f%15.10f\n", k, x[k], f[k], df[k]);
        return;
    }

    // Common tail (mode==1/3 entry falls through here).
    double dedx = zero;
    for (int k = 1; k <= minpq; ++k)
        dedx += theta[k - 1] * fock_interp[(k - 1) * norbs + (k - 1)];
    const double denow = -4.0 * dedx;
    const double enow = e;
    double x[13], f[13], df[13];
    x[ipoint] = xnow;
    f[ipoint] = enow;
    df[ipoint] = denow;
    double xmin = 0.0;
    spline(x, f, df, xhigh, xlow, xmin, npnts);
    if (eold_ref - enow <= ff * (eold_ref - emin) && ipoint <= 10) {
        xnow = xmin;
        for (int c = 0; c < norbs * norbs; ++c) cp[c] = vec_interp[c];
        for (int k = 1; k <= minpq; ++k) {
            const double ck = std::cos(xnow * theta[k - 1]);
            const double sk = std::sin(xnow * theta[k - 1]);
            if (debug)
                std::printf(" ROTATION ANGLE:%12.4f\n", sk * 57.29578);
            for (int r = 1; r <= norbs; ++r) {
                cp[(k - 1) * norbs + (r - 1)] =
                    ck * vec_interp[(k - 1) * norbs + (r - 1)] -
                    sk * vec_interp[(np + k - 1) * norbs + (r - 1)];
                cp[(np + k - 1) * norbs + (r - 1)] =
                    sk * vec_interp[(k - 1) * norbs + (r - 1)] +
                    ck * vec_interp[(np + k - 1) * norbs + (r - 1)];
            }
        }
        mode = 3;
        return;
    }
    if (mode != 2) {
        for (int c = 0; c < norbs * norbs; ++c) vec_interp[c] = cp[c];
        mode = 2;
    }
    eold_ref = enow;
    if (npnts <= 200) return;
    std::printf("  K     X(K)        F(K)         DF(K)\n");
    for (int k = 1; k <= npnts; ++k)
        std::printf("%3d%10.5f%15.10f%15.10f\n", k, x[k], f[k], df[k]);
}
