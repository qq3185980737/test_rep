// linmin.cpp — C++ translation of MOPAC 2016 "linmin.F90".
// Line minimisation: parabolic fit over three points (vt/phi), Thiel's
// alpha estimate, step limiting, and best-point restoration via exchng.
#include "linmin.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;

extern void compfg(const std::vector<double>& xparam, bool int_flag,
                   double& escf, bool fulscf, std::vector<double>& grad,
                   bool lgrad);
extern void exchng(double a, double& b, double c, double& d, double t,
                   double& q, const std::vector<double>& x,
                   std::vector<double>& y, int n);

void linmin(double* xparam, double& alpha, double* pvect, int nvar,
            double& funct, bool& okf, int& ic, double dott) {
    static bool print = false;
    static double alpold = 0.0, xmaxm = 0.0, delta1 = 0.0, delta2 = 0.0,
                  xnear = 0.0, ymaxst = 0.0;
    static int icalcn = 0, maxlin = 0;

    if (icalcn != numcal) {
        xmaxm = 0.4;
        delta2 = 0.001;
        delta1 = (keywrd.find("NOTH") == std::string::npos) ? 0.5 : 0.1;
        alpha = 1.0;
        maxlin = 15;
        xnear = 1e-4;
        if (nvar == 1) {
            pvect[0] = 0.01;
            alpha = 1.0;
            xnear = 1e-5;
            delta1 = 0.0005;
            delta2 = 0.0001;
            if (keywrd.find("PREC") != std::string::npos) {
                delta1 = 0.0;
                delta2 *= 0.01;
            }
            maxlin = 50;
        }
        cosine = 99.99;
        ymaxst = 0.4;
        print = keywrd.find("LINMIN") != std::string::npos;
        icalcn = numcal;
    }

    int nsame = 0;
    std::vector<double> xparef(nvar + 1), xstor(nvar + 1), xp(nvar + 1),
        grad(nvar + 1, 0.0);
    for (int i = 1; i <= nvar; ++i) xparef[i] = xparam[i - 1];

    xmaxm = 0.0;
    for (int i = 1; i <= nvar; ++i)
        xmaxm = std::max(xmaxm, std::abs(pvect[i - 1]));
    xmaxm = ymaxst / xmaxm;

    double fin = funct;
    double ssqlst = funct;
    bool diis = (ic == 1 && nvar > 1);
    double phi[4], vt[4];
    phi[1] = funct;
    alpha = 1.0;
    vt[1] = 0.0;
    vt[2] = alpha;
    vt[2] = std::min(xmaxm, vt[2]);
    double fmax = funct, fmin = funct;
    alpha = vt[2];
    for (int i = 1; i <= nvar; ++i) xp[i] = xparef[i] + alpha * pvect[i - 1];
    compfg(xp, true, phi[2], true, grad, false);
    fmax = std::max(phi[2], fmax);
    fmin = std::min(phi[2], fmin);
    double sqstor = 0.0, energy = 0.0, estor = 0.0, alfs = 0.0;
    exchng(phi[2], sqstor, energy, estor, alpha, alfs, xp, xstor, nvar);

    if (!diis) {
        if (nvar > 1) {
            alpha = -alpha * alpha * dott /
                    (2.0 * (phi[2] - ssqlst - alpha * dott));
            alpha = std::min(2.0, alpha);
        } else {
            if (phi[2] < phi[1]) alpha = 2.0 * alpha;
            else alpha = -alpha;
        }
        okf = okf || (phi[2] < ssqlst);
        if (delta1 > 0.3) {
            if (okf && alpha < 2.0) goto l190;
        }
        vt[3] = alpha;
        int left, center, right;
        if (vt[3] <= 1.0) {
            left = 3; center = 1; right = 2;
        } else {
            left = 1; center = 2; right = 3;
        }
        for (int i = 1; i <= nvar; ++i)
            xp[i] = xparef[i] + alpha * pvect[i - 1];
        compfg(xp, true, funct, true, grad, false);
        fmax = std::max(funct, fmax);
        fmin = std::min(funct, fmin);
        if (funct < sqstor)
            exchng(funct, sqstor, energy, estor, alpha, alfs, xp, xstor, nvar);
        okf = okf || (funct < fin);
        phi[3] = funct;
        if (print)
            std::printf(" ---QLINMN LEFT %.8f %.8f CENTER %.8f %.8f RIGHT %.8f %.8f\n",
                        vt[1], phi[1] - fin, vt[2], phi[2] - fin, vt[3], phi[3] - fin);
        alpold = 0.0;

        for (int ictr = 3; ictr <= maxlin; ++ictr) {
            double a = vt[2] - vt[3];
            double b = vt[3] - vt[1];
            double gamma = vt[1] - vt[2];
            if (std::abs(a * b * gamma) > xnear) {
                a = -(phi[1] * a + phi[2] * b + phi[3] * gamma) / (a * b * gamma);
            } else {
                break;  // two points very close together
            }
            b = (phi[1] - phi[2]) / gamma - a * (vt[1] + vt[2]);
            double s;
            if (a <= 0.0) {
                if (phi[right] <= phi[left])
                    a = 3.0 * vt[right] - 2.0 * vt[center];
                else
                    a = 3.0 * vt[left] - 2.0 * vt[center];
                s = a - alpold;
                if (std::abs(s) > xmaxm)
                    s = (xmaxm * (s >= 0 ? 1.0 : -1.0)) * (1.0 + 0.01 * (xmaxm / s));
                a = s + alpold;
            } else {
                a = -b / (2.0 * a);
                s = a - alpold;
                double xxm = 2.0 * xmaxm;
                if (std::abs(s) > xxm)
                    s = (xxm * (s >= 0 ? 1.0 : -1.0)) * (1.0 + 0.01 * (xxm / s));
                a = s + alpold;
            }
            bool near = false;
            for (int i = 1; i <= 3; ++i) {
                if (!(std::abs(a - vt[i]) < delta1 * (1.0 + vt[i]) && okf))
                    continue;
                near = true;
                break;
            }
            if (near) break;
            for (int i = 1; i <= nvar; ++i)
                xp[i] = xparef[i] + a * pvect[i - 1];
            double funold = funct;
            compfg(xp, true, funct, true, grad, false);
            fmax = std::max(funct, fmax);
            fmin = std::min(funct, fmin);
            if (funct < sqstor)
                exchng(funct, sqstor, energy, estor, a, alfs, xp, xstor, nvar);
            okf = okf || (funct < fin);
            if (print)
                std::printf(" LINMIN NEW %.8f %.8f %.8f\n", a, funct, funct - fin);
            if (std::abs(funold - funct) < delta2 && okf) break;
            alpold = a;
            bool to140 =
                (a > vt[right] ||
                 (a > vt[center] && funct < phi[center]) ||
                 (a > vt[left] && a < vt[center] && funct > phi[center]));
            if (to140 || nsame > 4) {
                nsame = 0;
                vt[left] = a;
                phi[left] = funct;
            } else {
                nsame = nsame + 1;
                vt[right] = a;
                phi[right] = funct;
            }
            if (vt[center] >= vt[right]) std::swap(center, right);
            if (vt[left] >= vt[center]) std::swap(left, center);
            if (vt[center] < vt[right]) continue;
            std::swap(center, right);
        }
    }
l190:
    ic = 2;
    if (std::abs(estor - energy) < 1e-12) ic = 1;
    double hlast = funct - fin;
    double drop = sqstor - fin;
    exchng(sqstor, funct, estor, energy, alfs, alpha, xstor, xp, nvar);
    if (hlast != 0.0 && hlast > 0.5 * drop)
        compfg(xp, true, funct, true, grad, false);
    okf = funct < ssqlst || diis;
    if (funct >= ssqlst) {
        for (int i = 1; i <= nvar; ++i) xparam[i - 1] = xp[i];
        return;
    }
    if (alpha < 0.0) {
        alpha = -alpha;
        for (int i = 1; i <= nvar; ++i) pvect[i - 1] = -pvect[i - 1];
    }
    for (int i = 1; i <= nvar; ++i) xparam[i - 1] = xp[i];
}
