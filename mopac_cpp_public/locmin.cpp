// locmin.cpp — C++ translation of MOPAC 2016 "locmin.F90".
// Line search for NLLSQ: minimises |gradient|^2 along direction p by
// parabolic fit over vt/phi triple with step limiting (xmaxm *= 3/iter)
// and early exits on tiny projected steps or false minima.
#include "locmin.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "molkst_C.h"

using namespace molkst_C;

extern void compfg(const std::vector<double>& xparam, bool int_flag,
                   double& escf, bool fulscf, std::vector<double>& grad,
                   bool lgrad);
extern void exchng(double a, double& b, double c, double& d, double t,
                   double& q, const std::vector<double>& x,
                   std::vector<double>& y, int n);
extern double ddot(int n, const double* x, int incx, const double* y,
                   int incy);

void locmin(int m, double* xparam, int nvar, double* p, double& ssq,
            double& alf, double* efs, int& ncount) {
    static double xmaxm = 0.0, scale = 0.0, eps = 0.0, tee = 0.0, ymaxst = 0.0,
                  xcrit = 0.0;
    static int mxcnt2 = 0, iprint = 0, icalcn = 0;
    static bool debug = false;

    if (icalcn != numcal) {
        icalcn = numcal;
        xmaxm = 1e9;
        scale = 1.0;
        eps = 1e-5;
        debug = keywrd.find("LINMIN") != std::string::npos;
        tee = 1e-2;
        ymaxst = 0.005;
        xcrit = 0.0002;
        mxcnt2 = 30;
        iprint = 0;
        if (debug) iprint = -1;
    }

    xmaxm = 1e-11;
    for (int i = 1; i <= nvar; ++i)
        xmaxm = std::max(xmaxm, std::abs(p[i - 1]));
    double xminm = xmaxm * scale;
    xmaxm = ymaxst / xmaxm / scale;

    std::vector<double> xp(nvar + 1), xstor(nvar + 1), gstor(m + 1, 0.0),
        efv(nvar + 1, 0.0), pv(nvar + 1, 0.0);
    for (int i = 1; i <= nvar; ++i) {
        xp[i] = xparam[i - 1];
        pv[i] = p[i - 1];
    }

    double fin = ssq;
    bool lower = false;
    double t = alf;
    double phi[4], vt[4];
    phi[1] = ssq;
    vt[1] = 0.0;
    vt[2] = t / 4.0;
    vt[2] = std::min(xmaxm, vt[2]);
    t = vt[2];
    double tscale = t * scale;
    for (int i = 1; i <= nvar; ++i) xp[i] += pv[i] * tscale;
    compfg(xp, true, escf, true, efv, true);
    phi[2] = ddot(nvar, &efv[1], 1, &efv[1], 1);
    double sqstor = 0.0, energy = 0.0, estor = 0.0, alfs = 0.0;
    exchng(phi[2], sqstor, energy, estor, t, alfs, xp, xstor, nvar);
    for (int i = 1; i <= m; ++i) gstor[i] = efv[i];

    int left, center, right;
    if (phi[1] <= phi[2]) {
        vt[3] = -vt[2];
        left = 3; center = 1; right = 2;
    } else {
        vt[3] = 2.0 * vt[2];
        left = 1; center = 2; right = 3;
    }
    double tlast = vt[3];
    t = tlast - t;
    tscale = t * scale;
    for (int i = 1; i <= nvar; ++i) xp[i] += pv[i] * tscale;
    double flast = phi[2];
    compfg(xp, true, escf, true, efv, true);
    double f = ddot(nvar, &efv[1], 1, &efv[1], 1);
    if (f < sqstor)
        exchng(f, sqstor, energy, estor, t, alfs, xp, xstor, nvar);
    for (int i = 1; i <= m; ++i) gstor[i] = efv[i];
    if (f < fin) lower = true;
    ncount += 2;
    phi[3] = f;
    if (iprint < 0)
        std::printf(" ---LOCMIN LEFT %.6f %.6f CENTER %.6f %.6f RIGHT %.6f %.6f\n",
                    vt[1], std::sqrt(phi[1]), vt[2], std::sqrt(phi[2]),
                    vt[3], std::sqrt(phi[3]));

    int mxct = mxcnt2;
    for (int ictr = 3; ictr <= mxct; ++ictr) {
        xmaxm = xmaxm * 3.0;
        double a = vt[2] - vt[3];
        double b = vt[3] - vt[1];
        double gamma = vt[1] - vt[2];
        if (a == 0.0) a = 1e-20;
        if (b == 0.0) b = 1e-20;
        if (gamma == 0.0) gamma = 1e-20;
        double abg = -(phi[1] * a + phi[2] * b + phi[3] * gamma) / a;
        abg = abg / b;
        abg = abg / gamma;
        a = abg;
        b = (phi[1] - phi[2]) / gamma - a * (vt[1] + vt[2]);
        double s;
        if (a <= 0.0) {
            if (phi[right] <= phi[left])
                t = 3.0 * vt[right] - 2.0 * vt[center];
            else
                t = 3.0 * vt[left] - 2.0 * vt[center];
            s = t - tlast;
            t = s + tlast;
        } else {
            t = -b / (2.0 * a);
            s = t - tlast;
            double amdis;
            if (s <= 0.0) {
                if (s == 0.0) break;
                amdis = vt[left] - tlast - xmaxm;
            } else {
                amdis = vt[right] - tlast + xmaxm;
            }
            if (std::abs(s) > std::abs(amdis)) s = amdis;
            t = s + tlast;
        }
        if (ictr > 3 && std::abs(s * xminm) < xcrit) break;
        t = s + tlast;
        double sscale = s * scale;
        for (int i = 1; i <= nvar; ++i) xp[i] += pv[i] * sscale;
        flast = f;
        compfg(xp, true, escf, true, efv, true);
        f = ddot(nvar, &efv[1], 1, &efv[1], 1);
        if (f < sqstor)
            exchng(f, sqstor, energy, estor, t, alfs, xp, xstor, nvar);
        for (int i = 1; i <= m; ++i) gstor[i] = efv[i];
        if (f < fin) lower = true;
        ncount += 1;
        if (iprint < 0)
            std::printf(" LOCMIN NEW %.6f %.6f %.6f %.6f %.6f\n",
                        vt[left], std::sqrt(phi[left]), vt[center],
                        std::sqrt(phi[center]), t, std::sqrt(f));
        // Test for excited states and potholes.
        if (std::abs(vt[center]) > 1e-10) goto l200;
        if (std::abs(t) / (std::abs(vt[left]) + 1e-15) > 0.3333) goto l200;
        if (2.5 * f - phi[right] - phi[left] < 0.5 * phi[center]) goto l200;
        break;  // stuck on a false minimum
    l200:;
        // Main stopping tests.
        if (debug) std::printf(" F/FLAST %.6f\n", f / flast);
        if (lower && f / flast > 0.995) {
            if (std::abs(t - tlast) <= eps * std::abs(t + tlast) + tee) break;
            double sum = std::min(std::min(std::abs(f - phi[1]),
                                           std::abs(f - phi[2])),
                                  std::abs(f - phi[3]));
            double sum2 = (fin - sqstor) * 0.05;
            if (sum < sum2) break;
        }
        tlast = t;
        bool store_right = !(t > vt[right] || (t > vt[center] && f < phi[center]) ||
                           (t > vt[left] && t < vt[center] && f > phi[center]));
        if (store_right) {
            vt[right] = t;
            phi[right] = f;
        } else {
            vt[left] = t;
            phi[left] = f;
        }
        if (vt[center] >= vt[right]) std::swap(center, right);
        if (vt[left] >= vt[center]) std::swap(left, center);
        if (vt[center] < vt[right]) continue;
        std::swap(center, right);
    }

    exchng(sqstor, f, estor, energy, alfs, t, xstor, xp, nvar);
    for (int i = 1; i <= m; ++i) efs[i - 1] = gstor[i];
    ssq = f;
    alf = t;
    if (t < 0.0) {
        t = -t;
        for (int i = 1; i <= nvar; ++i) p[i - 1] = -p[i - 1];
    }
    alf = t;
    for (int i = 1; i <= nvar; ++i) xparam[i - 1] = xp[i];
}
