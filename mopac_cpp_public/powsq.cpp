// powsq.cpp — C++ translation of "powsq.F90".
// POWSQ geometry optimizer: minimises the gradient norm (eigenvector-following-like).
// Includes internal line search "search" and restart I/O "powsav".

#include "powsq.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

// Forward declarations of external routines (ported elsewhere or stubbed).
extern "C" double ddot_(int* n, const double* dx, int* incx,
                        const double* dy, int* incy);
static inline double ddot(int n, const double* x, int, const double* y, int) {
    int incx = 1, incy = 1;
    return ddot_(&n, x, &incx, y, &incy);
}
double second(int);
double reada(const std::string& line, int pos);
void compfg(const std::vector<double>& xparam, bool, double& escf, bool,
            std::vector<double>& grad, bool);
void rsp(const std::vector<double>& a, int n, std::vector<double>& eig,
         std::vector<double>& pvec);
void vecprt(const std::vector<double>& a, int n);
void prttim(double tleft, double& tprt, char& txt);
void geout(int iw);
void den_in_out(int);
void mopend(const std::string& s);
void to_screen(const std::string& s);

// maps_C / ef_C / meci_C — minimal forward decls for powsav.
namespace maps_C { extern int latom; }
namespace ef_C {
extern double alparm[3];
extern double x0, x1, x2;
}
namespace meci_C { extern int jloop; }

using common_arrays_C::grad;
using common_arrays_C::gnext1;
using common_arrays_C::gmin1;
using common_arrays_C::loc;
using common_arrays_C::geo;
using common_arrays_C::xparam;
using common_arrays_C::aicorr;

using molkst_C::escf;
using molkst_C::keywrd;
using molkst_C::last;
using molkst_C::line;
using molkst_C::nscf;
using molkst_C::nvar;
using molkst_C::numat;
using molkst_C::norbs;
using molkst_C::tdump;
using molkst_C::tleft;
using molkst_C::time0;
using molkst_C::numcal;

namespace {

int icalcn = 0;
bool debug = false, restrt = false, times = false, scf1 = false, resfil = false;
double time1 = 0, time2 = 0, tlast = 0, xinc = 0, rho2 = 0, tol2 = 0;
int icyc = 0, iloop = 1;

void powsav(std::vector<double>& hess, std::vector<double>& gradv,
            std::vector<double>& xp, std::vector<double>& pmat,
            int& iloop_v, std::vector<double>& bmat, int* ipow);

// Line search: minimise gradient norm along direction sig.
void search_powsq(std::vector<double>& xp, double& alpha,
                  const std::vector<double>& sig, int nvar_,
                  double& gmin, double& funct, double& amin, double& anext) {
    static bool debug_s = false;
    static double g = 100.0, tiny = 0.1;
    static int looks = 0, icalcn_s = 0;
    static double tolerg = 0.02;

    if (icalcn_s != numcal) {
        icalcn_s = numcal;
        debug_s = keywrd.find("LINMIN") != std::string::npos;
        looks = 0; tiny = 0.1; tolerg = 0.02; g = 100.0;
        alpha = 0.1;
    }
    std::vector<double> gref(nvar_ + 1), xmin1(nvar_ + 1), xref(nvar_ + 1);
    std::vector<double> gradv(nvar_ + 1, 0.0);
    for (int i = 1; i <= nvar_; ++i) {
        gref[i] = gmin1[i];
        xmin1[i] = xp[i];
        xref[i] = xp[i];
    }
    for (int i = 1; i <= nvar_; ++i) gnext1[i] = gmin1[i];
    if (std::abs(alpha) > 0.2) alpha = std::copysign(0.2, alpha);

    double gb = ddot(nvar_, &gref[1], 1, &gmin1[1], 1);
    double gstore = gb;
    double amin_ = 0, gminn = 1e9, ta = 0, tb = 0;
    double ga = gb; gb = 1e9;
    int itrys = 0;
    // goto 30:
    auto step30 = [&]() {
        for (int i = 1; i <= nvar_; ++i) xp[i] = xref[i] + alpha * sig[i];
        std::fill(gradv.begin(), gradv.end(), 0.0);
        compfg(xp, true, funct, true, gradv, true);
        ++looks;
        g = ddot(nvar_, &gref[1], 1, &gradv[1], 1);
        double gtot = std::sqrt(ddot(nvar_, &gradv[1], 1, &gradv[1], 1));
        if (gtot < gminn) {
            gminn = gtot;
            if (std::abs(amin_ - alpha) > 1e-2) {
                anext = amin_;
                for (int i = 1; i <= nvar_; ++i) gnext1[i] = gmin1[i];
            }
            amin_ = alpha;
            if (gminn < gmin) {
                for (int i = 1; i <= nvar_; ++i) xmin1[i] = xp[i];
            }
            for (int i = 1; i <= nvar_; ++i) gmin1[i] = gradv[i];
            gmin = std::min(gminn, gmin);
        }
    };
    step30();  // first evaluation (itrys=0)
    while (true) {
        double sum = ga / (ga - gb);
        ++itrys;
        if (std::abs(sum) > 3.0) sum = std::copysign(3.0, sum);
        alpha = (tb - ta) * sum + ta;
        step30();
        if (itrys > 8) break;
        if (std::abs(g / gstore) < tiny || std::abs(g) < tolerg) break;
        if (std::abs(g) < std::max(std::abs(ga), std::abs(gb)) ||
            (ga * gb > 0.0 && g * ga < 0.0)) {
            if (std::abs(gb) < std::abs(ga)) { ta = alpha; ga = g; }
            else { tb = alpha; gb = g; }
        } else break;
    }
    gminn = std::sqrt(ddot(nvar_, &gmin1[1], 1, &gmin1[1], 1));
    for (int i = 1; i <= nvar_; ++i) xp[i] = xmin1[i];
    if (gminn > gmin)
        for (int i = 1; i <= nvar_; ++i) xp[i] = xref[i];
    amin = amin_;
}

}  // namespace

void powsq() {
    using chanel_C::iw;
    using chanel_C::log;

    tleft = tleft - second(1) + time0;
    double tstep = 0.0;
    if (icalcn != numcal) {
        icalcn = numcal;
        gnext1.assign(nvar + 1, 0.0);
        gmin1.assign(nvar + 1, 0.0);
        restrt = keywrd.find("RESTART") != std::string::npos;
        int maxcyc = 100000;
        size_t p = keywrd.find("CYCLES");
        if (p != std::string::npos) maxcyc = (int)reada(keywrd, (int)p);
        scf1 = keywrd.find("1SCF") != std::string::npos;
        time1 = second(2); time2 = time1;
        icyc = 0;
        times = keywrd.find("TIME") != std::string::npos;
        tlast = tleft;
        resfil = false;
        last = 0;
        iloop = 1;
        xinc = funcon_C::a0 * 0.01;
        rho2 = 1e-4;
        tol2 = 0.4;
        if (keywrd.find("PREC") != std::string::npos) tol2 = 1e-2;
        p = keywrd.find("GNORM");
        if (p != std::string::npos) {
            tol2 = reada(keywrd, (int)p);
            if (tol2 < 0.01 && keywrd.find("LET") == std::string::npos) tol2 = 0.01;
        }
        debug = keywrd.find("POWSQ") != std::string::npos;

        std::vector<double> hess((nvar + 1) * (nvar + 1), 0.0);
        std::vector<double> pmat(nvar * nvar + 1, 0.0);
        std::vector<double> bmat((nvar + 1) * (nvar + 1), 0.0);
        int ipow[10] = {};
        if (restrt) {
            ipow[9] = 0;
            powsav(hess, gmin1, xparam, pmat, iloop, bmat, ipow);
            icyc = ipow[3];
            if (scf1) {
                // goto 390
                std::fill(grad.begin(), grad.end(), 0.0);
                last = 1;
                compfg(xparam, true, escf, true, grad, true);
                for (int i = 1; i <= nvar; ++i) grad[i] = gmin1[i];
                molkst_C::iflepo = scf1 ? 13 : 11;
                return;
            }
            nscf = ipow[8];
            for (int i = 1; i <= nvar; ++i) { grad[i] = gmin1[i]; gnext1[i] = gmin1[i]; }
        }
    }
    nvar = std::abs(nvar);
    std::fill(grad.begin(), grad.end(), 0.0);
    compfg(xparam, true, escf, true, grad, true);
    double gmin = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
    for (int i = 1; i <= nvar; ++i) gnext1[i] = grad[i];
    for (int i = 1; i <= nvar; ++i) gmin1[i] = gnext1[i];

    // Hessian finite-difference loop.
    int nva = nvar;
    std::vector<double> hess((nva + 1) * (nva + 1), 0.0);
    std::vector<double> pmat(nva * nva + 1, 0.0);
    std::vector<double> bmat((nva + 1) * (nva + 1), 0.0);
    int ilpr = iloop;
    for (; iloop <= nvar; ++iloop) {
        time1 = second(1);
        xparam[iloop] += xinc;
        std::fill(grad.begin(), grad.end(), 0.0);
        compfg(xparam, true, escf, true, grad, true);
        if (scf1) {
            std::fill(grad.begin(), grad.end(), 0.0);
            last = 1;
            compfg(xparam, true, escf, true, grad, true);
            for (int i = 1; i <= nvar; ++i) grad[i] = gmin1[i];
            molkst_C::iflepo = scf1 ? 13 : 11;
            return;
        }
        grad[iloop] += 1e-5;
        xparam[iloop] -= xinc;
        for (int i = 1; i <= nvar; ++i)
            hess[iloop * (nva + 1) + i] = -(grad[i] - gnext1[i]) / xinc;
        if (chanel_C::iw0 >= 0) {
            line = std::to_string(iloop) + " of " + std::to_string(nvar) + " steps completed";
            to_screen(line);
        }
        time2 = second(2);
        tstep = time2 - time1;
        if (tlast - tleft > tdump) {
            tlast = tleft; resfil = true;
            int ipow[10] = {};
            ipow[9] = 2; ipow[3] = icyc; ipow[8] = nscf;
            int ii = iloop;
            powsav(hess, gmin1, xparam, pmat, ii, bmat, ipow);
        }
        if (tleft >= tstep * 2.0 && iloop - ilpr <= 100000) continue;
        int ipow[10] = {};
        ipow[9] = 1; ipow[8] = nscf; ipow[3] = icyc;
        int ii = iloop;
        powsav(hess, gmin1, xparam, pmat, ii, bmat, ipow);
        molkst_C::iflepo = -1;
        return;
    }

    // Scale Hessian.
    std::vector<double> work(nvar + 1);
    for (int i = 1; i <= nvar; ++i) {
        double sum = 0;
        for (int j = 1; j <= nvar; ++j) sum += hess[i * (nva + 1) + j] * hess[i * (nva + 1) + j];
        work[i] = 1.0 / std::sqrt(sum);
    }
    for (int i = 1; i <= nvar; ++i)
        for (int j = 1; j <= nvar; ++j) hess[i * (nva + 1) + j] *= work[i];
    for (int i = 1; i <= nvar; ++i) {
        for (int j = 1; j <= nvar; ++j) bmat[i * (nva + 1) + j] = 0.0;
        bmat[i * (nva + 1) + i] = work[i] * 2.0;
    }

    iloop = -99;
    tstep *= 4.0;
    int jcyc = icyc;

    // Main loop 130.
    std::vector<double> sig(nvar + 1), e1(nvar + 1), e2(nvar + 1),
        p(nvar + 1), eig(nvar + 1), q(nvar + 1);
    while (true) {
        if (tlast - tleft > tdump) {
            tlast = tleft; resfil = true;
            int ipow[10] = {};
            ipow[9] = 2; ipow[3] = icyc; ipow[8] = nscf;
            int ii = iloop;
            powsav(hess, gmin1, xparam, pmat, ii, bmat, ipow);
        }
        if (tleft < tstep * 2.0 || icyc - jcyc > 100000) {
            int ipow[10] = {};
            ipow[9] = 1; ipow[8] = nscf; ipow[3] = icyc;
            int ii = iloop;
            powsav(hess, gmin1, xparam, pmat, ii, bmat, ipow);
            molkst_C::iflepo = -1;
            return;
        }
        // 140.
        int ij = 0;
        for (int j = 1; j <= nvar; ++j)
            for (int i = 1; i <= j; ++i) {
                ++ij;
                double sum = 0;
                for (int k = 1; k <= nvar; ++k)
                    sum += hess[i * (nva + 1) + k] * hess[j * (nva + 1) + k];
                pmat[ij] = sum;
            }
        for (int i = 1; i <= nvar; ++i) {
            double sum = 0;
            for (int k = 1; k <= nvar; ++k) sum -= hess[i * (nva + 1) + k] * gmin1[k];
            p[i] = -sum;
        }
        std::vector<double> pvec(nvar * nvar + 1, 0.0);
        rsp(pmat, nvar, eig, pvec);
        int l = 0;
        if (eig[1] >= rho2) {
            ij = 0;
            for (int i = 1; i <= nvar; ++i)
                for (int j = 1; j <= i; ++j) {
                    ++ij;
                    double sum = 0;
                    for (int k = 1; k <= nvar; ++k)
                        sum += pvec[(k - 1) * nvar + j] * pvec[(k - 1) * nvar + i] / eig[k];
                    pmat[ij] = sum;
                }
            for (int i = 1; i <= nvar; ++i) {
                double sum = 0;
                for (int k = 1; k <= i; ++k) {
                    int ik = i * (i - 1) / 2 + k;
                    sum += pmat[ik] * p[k];
                }
                for (int k = i + 1; k <= nvar; ++k) {
                    int ik = k * (k - 1) / 2 + i;
                    sum += pmat[ik] * p[k];
                }
                q[i] = sum;
            }
        } else {
            for (int i = 1; i <= nvar; ++i) q[i] = pvec[i];
        }
        for (int i = 1; i <= nvar; ++i) {
            sig[i] = 0;
            for (int j = 1; j <= nvar; ++j) sig[i] += q[j] * bmat[i * (nva + 1) + j];
        }
        double alpha = 0, funct = 0, amin = 0, anext = 0;
        search_powsq(xparam, alpha, sig, nvar, gmin, funct, amin, anext);
        if (nvar == 1) {
            std::fill(grad.begin(), grad.end(), 0.0);
            last = 1;
            compfg(xparam, true, escf, true, grad, true);
            for (int i = 1; i <= nvar; ++i) grad[i] = gmin1[i];
            molkst_C::iflepo = scf1 ? 13 : 11;
            return;
        }
        double rmx = 0;
        for (int k = 1; k <= nvar; ++k) rmx = std::max(std::abs(gmin1[k]), rmx);
        if (rmx < tol2) {
            std::fill(grad.begin(), grad.end(), 0.0);
            last = 1;
            compfg(xparam, true, escf, true, grad, true);
            for (int i = 1; i <= nvar; ++i) grad[i] = gmin1[i];
            molkst_C::iflepo = scf1 ? 13 : 11;
            return;
        }
        for (int i = 1; i <= nvar; ++i)
            e1[i] = (gmin1[i] - gnext1[i]) / (amin - anext);
        double rmu = ddot(nvar, &e1[1], 1, &gmin1[1], 1) /
                     ddot(nvar, &gmin1[1], 1, &gmin1[1], 1);
        for (int i = 1; i <= nvar; ++i) e2[i] = e1[i] - rmu * gmin1[i];
        double sk = 1.0 / std::sqrt(ddot(nvar, &e2[1], 1, &e2[1], 1));
        for (int i = 1; i <= nvar; ++i) { sig[i] *= sk; e2[i] *= sk; }
        double pmax = -1e20; int id_ = 1;
        for (int i = 1; i <= nvar; ++i) {
            if (std::abs(p[i] * q[i]) <= pmax) continue;
            pmax = std::abs(p[i] * q[i]);
            id_ = i;
        }
        for (int j = 1; j <= nvar; ++j) hess[id_ * (nva + 1) + j] = -e2[j];
        for (int i = 1; i <= nvar; ++i) bmat[i * (nva + 1) + id_] = sig[i] / funcon_C::a0;
        for (int i = 1; i <= nvar; ++i) gnext1[i] = gmin1[i];
        time1 = time2;
        time2 = second(2);
        tstep = time2 - time1;
        tleft -= tstep;
        if (tleft < 0) tleft = -0.1;
        ++icyc;
        double tprt = 0; char txt = ' ';
        prttim(tleft, tprt, txt);
        if (resfil) {
            line = "  RESTART FILE WRITTEN";
            resfil = false;
        } else {
            line = " CYCLE: " + std::to_string(icyc);
        }
        to_screen(line);
    }
}

namespace {
void powsav(std::vector<double>& hess, std::vector<double>& gradv,
            std::vector<double>& xp, std::vector<double>& pmat,
            int& iloop_v, std::vector<double>& bmat, int* ipow) {
    // Restart I/O: in C++ this uses a simple binary stream.
    std::fstream fs(chanel_C::restart_fn,
                    std::ios::binary | std::ios::in | std::ios::out);
    if (!fs) fs.open(chanel_C::restart_fn, std::ios::binary | std::ios::out);
    if (ipow[9] != 0) {
        for (int i = 1; i <= nvar; ++i) {
            int k = loc[1][i], l = loc[2][i];
            geo[l][k] = xp[i];
        }
        // Write record.
        fs.write(reinterpret_cast<const char*>(&numat), sizeof(numat));
        fs.write(reinterpret_cast<const char*>(&norbs), sizeof(norbs));
        fs.write(reinterpret_cast<const char*>(&xp[1]), nvar * sizeof(double));
        fs.write(reinterpret_cast<const char*>(ipow + 1), 9 * sizeof(int));
        fs.write(reinterpret_cast<const char*>(&iloop_v), sizeof(iloop_v));
        fs.write(reinterpret_cast<const char*>(&gradv[1]), nvar * sizeof(double));
        fs.write(reinterpret_cast<const char*>(&hess[1]), (nvar + 1) * (nvar + 1) * sizeof(double));
        fs.write(reinterpret_cast<const char*>(&bmat[1]), (nvar + 1) * (nvar + 1) * sizeof(double));
        int linear = nvar * (nvar + 1) / 2;
        fs.write(reinterpret_cast<const char*>(&pmat[1]), linear * sizeof(double));
    } else {
        int old_numat, old_norbs;
        fs.read(reinterpret_cast<char*>(&old_numat), sizeof(old_numat));
        fs.read(reinterpret_cast<char*>(&old_norbs), sizeof(old_norbs));
        fs.read(reinterpret_cast<char*>(&xp[1]), nvar * sizeof(double));
        if (norbs != old_norbs || numat != old_numat) {
            mopend("Restart file read in does not match current data set");
            return;
        }
        fs.read(reinterpret_cast<char*>(ipow + 1), 9 * sizeof(int));
        fs.read(reinterpret_cast<char*>(&iloop_v), sizeof(iloop_v));
        fs.read(reinterpret_cast<char*>(&gradv[1]), nvar * sizeof(double));
        fs.read(reinterpret_cast<char*>(&hess[1]), (nvar + 1) * (nvar + 1) * sizeof(double));
        fs.read(reinterpret_cast<char*>(&bmat[1]), (nvar + 1) * (nvar + 1) * sizeof(double));
        int linear = nvar * (nvar + 1) / 2;
        fs.read(reinterpret_cast<char*>(&pmat[1]), linear * sizeof(double));
        ++iloop_v;
    }
}
}  // namespace
