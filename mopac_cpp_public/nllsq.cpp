// nllsq.cpp — C++ translation of MOPAC 2016 "nllsq.F90".
// Non-derivative nonlinear least-squares minimizer: Levenberg-Marquardt
// with QR decomposition, Broyden rank-1 update, and orthogonal sidesteps.

#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "mopend.h"

// External dependency routines (ported elsewhere or stubbed).
extern double second(int);
double reada(const std::string& s, int pos);
void compfg(const std::vector<double>& xparam, bool, double& escf, bool,
              std::vector<double>& grad, bool);
void locmin(int m, double* xparam, int nvar, double* p, double& ssq,
                  double& alf, double* efs, int& ncount);
void parsav(int mode, int& nvar, int& m, double* q, double* r, double* grad,
            double* xlast, int* iiium);
void geout(int);
void prttim(double tleft, double& tprt, char& txt);
void to_screen(const std::string& line);
double ddot(int n, const double* dx, int incx, const double* dy, int incy);

using namespace common_arrays_C;
using namespace molkst_C;

namespace {
double sgn(double x, double y) { return y >= 0 ? std::fabs(x) : -std::fabs(x); }
}

void nllsq() {
    int m = nvar;
    std::vector<std::vector<double>> q(m + 1, std::vector<double>(m + 1, 0.0));
    std::vector<std::vector<double>> r(m + 3, std::vector<double>(m + 3, 0.0));
    std::vector<int> iiium(7, 0);
    std::vector<double> y(nvar + 1), efs(nvar + 1), p(nvar + 1), xlast(nvar + 1);

    bool middle = keywrd.find(" RESTART") != std::string::npos;
    int maxcyc = 100000;
    {
        size_t pos = keywrd.find(" CYCLES");
        if (pos != std::string::npos) maxcyc = (int)std::lround(reada(keywrd, (int)pos));
    }
    iflepo = 10;
    double pn = 0.0, ssq = 0.0, alf = 0.0;

    double tol2 = 0.4;
    {
        size_t pos = keywrd.find(" GNORM");
        if (pos != std::string::npos) {
            tol2 = reada(keywrd, (int)pos);
            if (tol2 < 0.01 && keywrd.find(" LET") == std::string::npos) tol2 = 0.01;
        }
    }
    last = 0;
    double tols1 = 1e-12, tols2 = 1e-10, tols5 = 1e-6, tols6 = 1e-3;
    int nrst = 4;
    double tlast = tleft;
    bool resfil = false;
    tleft = tleft - second(2) + time0;

    int ifrtl = 0, nsst = 0;
    int ixso = nvar, np1 = nvar + 1, np2 = nvar + 2;
    int icyc = 0, jcyc = 0, irst = 0, jrst = 1;
    double eps = tols5, t = tols6;
    for (int i = 1; i <= nvar; ++i) grad[i] = 0.0;

    int ncount = 1;
    if (middle) {
        for (int i = 1; i <= m; ++i) for (int j = 1; j <= m; ++j) q[i][j] = 0;
        for (int i = 1; i <= np2; ++i) for (int j = 1; j <= np2; ++j) r[i][j] = 0;
        parsav(0, nvar, m, &q[0][0], &r[0][0], grad.data(), xlast.data(), iiium.data());
        if (moperr) return;
        nscf = iiium[1];
        ncount = iiium[5];
        for (int i = 1; i <= nvar; ++i) xparam[i] = xlast[i];
        double time1 = second(2); (void)time1;
        if (keywrd.find(" 1SCF") != std::string::npos) { iflepo = 13; last = 1; return; }
        jcyc = icyc;
        goto main_loop;
    }

    compfg(xparam, true, escf, true, grad, true);
    ssq = ddot(nvar, grad.data(), 1, grad.data(), 1);
    ncount = 1;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= nvar; ++j) { r[i][j] = (i == j) ? 1.0 : 0.0; }
        for (int j = i; j <= m; ++j) { q[i][j] = 0; q[j][i] = 0; if (i == j) q[i][i] = 1.0; }
    }
    double temp = 0.0;
    for (int i = 1; i <= nvar; ++i) temp += xparam[i] * xparam[i];
    alf = 100.0 * (eps * std::sqrt(temp) + t);

main_loop:
    {
        double time1 = second(2); (void)time1;
        ++ifrtl; ++icyc; ++irst;
        double prt = std::sqrt(ssq);
        if (irst < nrst) {
            nsst = 0;
            for (int i = 1; i <= m; ++i) {
                temp = 0.0;
                for (int j = 1; j <= m; ++j) temp += q[j][i] * grad[j];
                efs[i] = -temp;
            }
            for (int j = 1; j <= nvar; ++j) {
                int jj = np1 - j;
                for (int i = 1; i <= j; ++i) {
                    int ii = np2 - i;
                    r[ii][jj] = r[i][j];
                }
            }
            for (int i = 1; i <= nvar; ++i) {
                int i1 = i + 1;
                y[i] = prt;
                double efsss = 0.0;
                if (i < nvar) for (int z = i1; z <= nvar; ++z) y[z] = 0.0;
                for (int j = i; j <= nvar; ++j) {
                    int ii = np2 - j, jj = np1 - j;
                    if (std::fabs(y[j]) >= std::fabs(r[ii][jj]))
                        temp = y[j] * std::sqrt(1.0 + (r[ii][jj] / y[j]) * (r[ii][jj] / y[j]));
                    else
                        temp = r[ii][jj] * std::sqrt(1.0 + (y[j] / r[ii][jj]) * (y[j] / r[ii][jj]));
                    double sin = r[ii][jj] / temp, cos = y[j] / temp;
                    r[ii][jj] = temp;
                    temp = efs[j];
                    efs[j] = sin * temp + cos * efsss;
                    efsss = sin * efsss - cos * temp;
                    if (j >= nvar) break;
                    int j1 = j + 1;
                    for (int k = j1; k <= nvar; ++k) {
                        int jk = np1 - k;
                        temp = r[ii][jk];
                        r[ii][jk] = sin * temp + cos * y[k];
                        y[k] = sin * y[k] - cos * temp;
                    }
                }
            }
            p[nvar] = efs[nvar] / r[2][1];
            for (int i = nvar - 1; i >= 1; --i) {
                temp = efs[i];
                int k = i + 1;
                int ii = np2 - i;
                for (int z = k; z <= nvar; ++z) {
                    int jj = np1 - z;
                    temp -= r[ii][jj] * p[z];
                }
                int jj = np1 - i;
                p[i] = temp / r[ii][jj];
            }
        } else {
            ++jrst; ++nsst;
            if (nsst >= ixso) goto sidestep_fail;
            if (jrst > nvar) jrst = 2;
            irst = 0;
            double work = pn * (std::fabs(p[1]) + pn);
            double temp2 = p[jrst];
            p[1] = temp2 * (p[1] + sgn(pn, p[1]));
            for (int z = 2; z <= nvar; ++z) p[z] = temp2 * p[z];
            p[jrst] -= work;
        }

        double pnlast = pn;
        pn = 0.0; double pn2 = 0.0;
        for (int i = 1; i <= nvar; ++i) { pn += std::fabs(p[i]); pn2 += p[i] * p[i]; }
        if (pn < 1e-20) {
            geout(1);
            mopend("SYSTEM DOES NOT APPEAR TO BE OPTIMIZABLE.");
            goto cleanup;
        }
        pn2 = std::max(1e-20, pn2);
        pn = std::sqrt(pn2);
        alf = std::min(1e20, alf);
        if (icyc > 1) {
            alf = alf * 1e-20 * pnlast / pn;
            alf = std::min(1e10, alf);
            alf *= 1e20;
        }
        double ttmp = alf * pn;
        if (ttmp < 1e-4) alf = 0.001 / pn;

        for (int i = 1; i <= nvar; ++i) efs[i] = xparam[i];
        double ssqlst = ssq;
        for (int i = 1; i <= nvar; ++i) { efs[i] = 0.0; xlast[i] = xparam[i]; }
        locmin(m, xparam.data(), nvar, p.data(), ssq, alf, efs.data(), ncount);
        if (ssqlst < ssq) {
            for (int i = 1; i <= nvar; ++i) xparam[i] = xlast[i];
            irst = nrst;
            pn = pnlast;
            goto cycle_report;
        }

        for (int i = 1; i <= nvar; ++i) {
            temp = 0.0;
            for (int j = i; j <= nvar; ++j) temp += r[i][j] * p[j];
            y[i] = temp;
        }
        double work = alf * pn2;
        double yn = 0.0;
        for (int i = 1; i <= m; ++i) {
            temp = 0.0;
            for (int j = 1; j <= nvar; ++j) temp += q[i][j] * y[j];
            temp = efs[i] - grad[i] - alf * temp;
            grad[i] = efs[i];
            yn += temp * temp;
            efs[i] = temp / work;
        }
        yn = std::sqrt(yn) / work;

        for (int i = 1; i <= m; ++i) {
            temp = 0.0;
            for (int j = 1; j <= m; ++j) temp += q[j][i] * efs[j];
            y[i] = temp;
        }

        if (m > np1) {
            double cnst = 1e-12;
            for (int i = np1; i <= m; ++i) cnst = std::max(std::fabs(y[np1]), cnst);
            double ytail = 0.0;
            for (int i = np1; i <= m; ++i) { double d = y[i] / cnst; ytail += d * d; }
            ytail = std::sqrt(ytail) * cnst;
            double bet = (1e25 / ytail) / (ytail + std::fabs(y[np1]));
            y[np1] = sgn(ytail + std::fabs(y[np1]), y[np1]);
            for (int i = 1; i <= m; ++i) {
                double tmp = 0.0;
                for (int z = np1; z <= m; ++z) tmp += q[i][z] * y[z] * 1e-25;
                tmp = bet * tmp;
                for (int z = np1; z <= m; ++z) q[i][z] -= tmp * y[z];
            }
            y[np1] = ytail;
            m = np1;
        } else {
            m = nvar;
        }
        {
            int ii = m;
            int jj = ii + 1;
            while (ii > 0) {
                jj = ii + 1;
                if (y[jj] != 0.0) {
                    double sin, cos;
                    if (std::fabs(y[ii]) >= std::fabs(y[jj]))
                        temp = std::fabs(y[ii]) * std::sqrt(1.0 + (y[jj] / y[ii]) * (y[jj] / y[ii]));
                    else
                        temp = std::fabs(y[jj]) * std::sqrt(1.0 + (y[ii] / y[jj]) * (y[ii] / y[jj]));
                    cos = y[ii] / temp; sin = y[jj] / temp;
                    y[ii] = temp;
                    for (int k = 1; k <= m; ++k) {
                        double t2 = cos * q[k][ii] + sin * q[k][jj];
                        double w2 = -sin * q[k][ii] + cos * q[k][jj];
                        q[k][ii] = t2; q[k][jj] = w2;
                    }
                    if (ii <= nvar) {
                        r[jj][ii] = -sin * r[ii][ii];
                        r[ii][ii] = cos * r[ii][ii];
                        if (jj <= nvar) {
                            for (int k = jj; k <= nvar; ++k) {
                                double t2 = cos * r[ii][k] + sin * r[jj][k];
                                double w2 = -sin * r[ii][k] + cos * r[jj][k];
                                r[ii][k] = t2; r[jj][k] = w2;
                            }
                        }
                    }
                }
                --ii;
                if (ii <= 0) break;
            }
        }
        for (int z = 1; z <= nvar; ++z) r[1][z] += yn * p[z];
        int jend = np1;
        if (m == nvar) jend = nvar;
        for (int j = 2; j <= jend; ++j) {
            int i = j - 1;
            if (r[j][i] == 0.0) continue;
            if (std::fabs(r[i][i]) >= std::fabs(r[j][i]))
                temp = std::fabs(r[i][i]) * std::sqrt(1.0 + (r[j][i] / r[i][i]) * (r[j][i] / r[i][i]));
            else
                temp = std::fabs(r[j][i]) * std::sqrt(1.0 + (r[i][i] / r[j][i]) * (r[i][i] / r[j][i]));
            double cos = r[i][i] / temp, sin = r[j][i] / temp;
            r[i][i] = temp;
            if (j <= nvar) {
                for (int k = j; k <= nvar; ++k) {
                    double t2 = cos * r[i][k] + sin * r[j][k];
                    double w2 = -sin * r[i][k] + cos * r[j][k];
                    r[i][k] = t2; r[j][k] = w2;
                }
            }
            for (int k = 1; k <= m; ++k) {
                double t2 = cos * q[k][i] + sin * q[k][j];
                double w2 = -sin * q[k][i] + cos * q[k][j];
                q[k][i] = t2; q[k][j] = w2;
            }
        }

        temp = 0.0;
        for (int i = 1; i <= nvar; ++i) temp += xparam[i] * xparam[i];
        double tolx = tols1 * std::sqrt(temp) + tols2;
        if (std::sqrt(alf * pn2) <= tolx) { std::printf("TEST ON xparam SATISFIED\n"); goto finish; }
        if (ssq >= 2.0 * nvar) goto cycle_report;
        bool conv = true;
        for (int i = 1; i <= nvar; ++i) if (std::fabs(grad[i]) >= tol2) { conv = false; break; }
        if (conv) { std::printf("TEST ON SSQ SATISFIED\n"); goto finish; }

    cycle_report:
        {
            double time2 = second(2); (void)time2;
            // (time accounting abbreviated)
        }
        double tprt = 0; char txt = ' ';
        prttim(tleft, tprt, txt);
        std::string l;
        char buf[128];
        std::snprintf(buf, sizeof(buf), " CYCLE:%6d TIME GRAD.:%10.3f HEAT:%g\n",
                      icyc, std::sqrt(ssq), escf);
        l = buf;
        to_screen(l);
        if (tleft > 0) resfil = false;
        if (tlast - tleft > tdump) {
            tlast = tleft;
            resfil = true;
            for (int i = 1; i <= nvar; ++i) xlast[i] = xparam[i];
            iiium[1] = nscf;
            parsav(2, nvar, m, &q[0][0], &r[0][0], grad.data(), xlast.data(), iiium.data());
            if (moperr) goto cleanup;
        }
        if (tleft > 0 && icyc - jcyc < maxcyc) goto main_loop;
        iiium[5] = ncount;
        for (int i = 1; i <= nvar; ++i) xlast[i] = xparam[i];
        iiium[1] = nscf;
        parsav(1, nvar, m, &q[0][0], &r[0][0], grad.data(), xlast.data(), iiium.data());
        if (moperr) goto cleanup;
        iflepo = -1;
        goto cleanup;
    }

sidestep_fail:
    std::printf("ATTEMPT TO GO DOWNHILL UNSUCCESSFUL AFTER %d ORTHOGONAL SEARCHES\n", ixso);
finish:
    last = 1;
cleanup:
    return;
}