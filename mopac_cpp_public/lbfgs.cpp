// lbfgs.cpp — C++ translation of MOPAC 2016 "lbfgs.F90" (3203 lines).
//
// Full L-BFGS-B limited-memory quasi-Newton optimizer, faithful port:
//   lbfgs        — driver used by geometry optimization (run_mopac dispatch)
//   lbfsav       — restart-file save/restore
//   setulb       — L-BFGS-B entry wrapper (wa workspace partitioning)
//   mainlb       — L-BFGS-B main kernel (state machine over task)
//   active, bmv, cauchy, cmprlb, errclb, formk, formt, freev, hpsolb,
//   lnsrlb, matupd, prn1lb, prn2lb, prn3lb, projgr, subsm, dcsrch, dcstep,
//   dpmeps, dpofa, dtrsl
//
// Adaptation notes:
//  - 1-based Fortran indexing kept throughout: vectors are n+1, matrices are
//    column-major views (Mat2D) so BLAS calls on columns stay contiguous.
//  - F90 write(iw,...) goes to stdout; endfile/backspace are no-ops.
//  - F90 open/write of restart file (lbfsav) uses restart_fn on ires.

#include "lbfgs.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "compfg.h"
#include "den_in_out.h"
#include "ef_C.h"
#include "geo_diff.h"
#include "geochk.h"
#include "geout.h"
#include "memory_error.h"
#include "molkst_C.h"
#include "mopend.h"
#include "prtgra.h"
#include "prttim.h"
#include "reada.h"
#include "second.h"
#include "symtry.h"
#include "timer.h"
#include "to_screen.h"
#include "write_cell.h"

using namespace chanel_C;
using namespace common_arrays_C;
using namespace ef_C;
using namespace molkst_C;

namespace lbfgs_impl {

// Column-major 1-based matrix view over a contiguous block (F90 style).
struct Mat2D {
  std::vector<double> data;  // owned storage (default / rc ctor)
  double* ptr;               // view pointer into external storage (wa)
  int rows, cols;
  bool view;
  Mat2D() : rows(0), cols(0), ptr(nullptr), view(false) {}
  Mat2D(int r, int c) : data((size_t)r * c, 0.0), rows(r), cols(c), view(false) {
    ptr = data.empty() ? nullptr : &data[0];
  }
  // External-view constructor: element (i,j) at p[(j-1)*rows + (i-1)].
  Mat2D(double* p, int r, int c) : ptr(p), rows(r), cols(c), view(true) {}
  double* base() { return view ? ptr : &data[0]; }
  const double* base() const { return view ? ptr : &data[0]; }
  double& at(int i, int j) { return base()[(j - 1) * rows + (i - 1)]; }
  const double& at(int i, int j) const { return base()[(j - 1) * rows + (i - 1)]; }
  double* col(int j) { return base() + (j - 1) * rows; }
  const double* col(int j) const { return base() + (j - 1) * rows; }
};

// 1-based pointer BLAS helpers: element i is x[1 + (i-1)*inc].
inline double ddot1(int n, const double* x, int incx, const double* y, int incy) {
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += x[(size_t)i * incx] * y[(size_t)i * incy];
  return s;
}
inline void dcopy1(int n, const double* x, int incx, double* y, int incy) {
  for (int i = 0; i < n; ++i) y[(size_t)i * incy] = x[(size_t)i * incx];
}
inline void daxpy1(int n, double a, const double* x, int incx, double* y, int incy) {
  for (int i = 0; i < n; ++i) y[(size_t)i * incy] += a * x[(size_t)i * incx];
}
inline void dscal1(int n, double a, double* x, int incx) {
  for (int i = 0; i < n; ++i) x[(size_t)i * incx] *= a;
}

}  // namespace lbfgs_impl

using lbfgs_impl::Mat2D;
using lbfgs_impl::ddot1;
using lbfgs_impl::dcopy1;
using lbfgs_impl::daxpy1;
using lbfgs_impl::dscal1;

// ---- forward declarations (internal L-BFGS-B library) ----
void lbfsav(double tt0, int mode, std::vector<double>& wa, int nwa,
                   std::vector<int>& iwa, int niwa, std::string& task,
                   std::string& csave, std::vector<int>& lsave, std::vector<int>& isave,
                   std::vector<double>& dsave, int& nstep, double& escf);
void setulb(int n, int m, double* x, double* l, double* u, int* nbd,
                   double& f, double* g, double factr, double pgtol,
                   std::vector<double>& wa, std::vector<int>& iwa, std::string& task,
                   int iprint, std::string& csave, std::vector<int>& lsave,
                   std::vector<int>& isave, std::vector<double>& dsave);
static void mainlb(int n, int m, double* x, double* l, double* u, int* nbd,
                   double& f, double* g, double factr, double pgtol,
                   Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& ss, Mat2D& wt,
                   Mat2D& wn, Mat2D& snd, double* z, double* r, double* d,
                   double* t, double* wa, int* Index, int* iwhere, int* indx2,
                   std::string& task, int iprint, std::string& csave,
                   std::vector<int>& lsave, int* isave, double* dsave);
static void active(int n, double* l, double* u, int* nbd, double* x,
                   int* iwhere, int iprint, bool& prjctd, bool& cnstnd,
                   bool& boxed);
static void bmv(int m, Mat2D& sy, Mat2D& wt, int col, std::vector<double>& v,
                std::vector<double>& p, int& info);
static void cauchy(int n, double* x, double* l, double* u, int* nbd,
                   const double* g, std::vector<int>& iorder, int* iwhere,
                   double* t, double* d, double* xcp, int m, Mat2D& wy, Mat2D& ws,
                   Mat2D& sy, Mat2D& wt, double& theta, int col, int head,
                   double* p, double* c, double* wbp, double* v, int& nint,
                   int iprint, double sbgnrm, int& info, double epsmch);
static void cmprlb(int n, int m, const double* x, const double* g,
                   Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& wt, const double* z,
                   double* r, double* wa, const int* Index, double theta, int col,
                   int head, int nfree, bool cnstnd, int& info);
static void errclb(int n, int m, double factr, double* l, double* u, int* nbd,
                   std::string& task, int& info, int& k);
static void formk(int n, int nsub, const int* ind, int nenter, int ileave,
                  const int* indx2, int iupdat, bool updatd, Mat2D& wn,
                  Mat2D& wn1, int m, Mat2D& ws, Mat2D& wy, Mat2D& sy, double theta,
                  int col, int head, int& info);
static void formt(int m, Mat2D& wt, Mat2D& sy, Mat2D& ss, int col, double theta, int& info);
static void freev(int n, int& nfree, int* Index, int& nenter, int& ileave,
                  int* indx2, const int* iwhere, bool& wrk,
                  bool updatd, bool cnstnd, int iprint, int iter);
static void hpsolb(int n, double* t, std::vector<int>& iorder, int iheap);
static void lnsrlb(int n, double* l, double* u, int* nbd, double* x, double f,
                   double& fold, double& gd, double& gdold, double* g,
                   double* d, double* r, double* t, double* z, double& stp,
                   double& dnorm, double& dtd, double& xstep, double& stpmx,
                   int iter, int& ifun, int& iback, int& nfgv, int& info,
                   std::string& task, bool boxed, bool cnstnd, std::string& csave,
                   int* isave, double* dsave);
static void matupd(int n, int m, Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& ss,
                   double* d, double* r, int& itail,
                   int iupdat, int& col, int& head, double& theta, double rr,
                   double dr, double stp, double dtd);
static void prn1lb(int n, int m, double* l, double* u, double* x, int iprint,
                   int itfile, double epsmch);
static void prn2lb(int n, double* x, double f, double* g, int iprint, int itfile,
                   int iter, int nfgv, int nact, double sbgnrm, int nint,
                   std::string& word, int iword, int iback, double stp, double xstep);
static void prn3lb(int n, double* x, double f, const std::string& task, int iprint,
                   int info, int itfile, int iter, int nfgv, int nintol, int nskip,
                   int nact, double sbgnrm, double time, int nint,
                   const std::string& word, int iback, double stp, double xstep, int k,
                   double cachyt, double sbtime, double lnscht);
static void projgr(int n, double* l, double* u, int* nbd, double* x, double* g,
                   double& sbgnrm);
static void subsm(int n, int m, int nsub, const int* ind, double* l, double* u,
                  int* nbd, double* x, double* d, Mat2D& ws, Mat2D& wy,
                  double theta, int col, int head, int& iword, std::vector<double>& wv,
                  Mat2D& wn, int iprint, int& info);
static void dcsrch(double f, double g, double& stp, double ftol, double gtol,
                   double xtol, double stpmin, double stpmax, std::string& task,
                   int* isave, double* dsave);
static void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy,
                   double& dy, double& stp, double fp, double dp, bool& brackt,
                   double stpmin, double stpmax);
static double dpmeps();
static void dpofa(Mat2D& a, int lda, int n, int& info, int r0 = 1, int c0 = 1);
static void dtrsl(Mat2D& t, int ldt, int n, double* b, int job, int& info);

// ============================================================================
// lbfgs — driver (run_mopac dispatch target "LBFGS")
// ============================================================================
void lbfgs(std::vector<double>& xparam_arr, double& escf) {
    lbfgs(xparam_arr.data(), escf);
}
void lbfgs(double* xparam_arr, double& escf) {
  std::vector<double> xparam((size_t)nvar + 1, 0.0);
  for (int i = 1; i <= nvar; ++i) xparam[i] = xparam_arr[i];
  std::vector<double> best_xparam((size_t)nvar + 1, 0.0);
  std::vector<double> best_gradients((size_t)nvar + 1, 0.0);
  std::vector<double> bot((size_t)nvar + 1, 0.0), gold((size_t)nvar + 1, 0.0);
  std::vector<double> top((size_t)nvar + 1, 0.0), xold((size_t)nvar + 1, 0.0);
  std::vector<int> best_nc((size_t)natoms + 1, 0), nbd((size_t)nvar + 1, 0);
  double best_funct = 1.e10, best_gnorm = 1.e10;
  bool resfil = false;
  double tlast = 0.0;
  int nflush = 1;
  int m = 12;
  std::fill(grad.begin(), grad.end(), 0.0);
  int niwa = 3 * nvar;
  int nwa = 2 * nvar * m + 4 * nvar + 11 * m * m + 8 * m;
  std::vector<double> wa((size_t)nwa + 1, 0.0);
  std::vector<int> iwa((size_t)niwa + 1, 0);
  std::string task = " Unused", csave = " Unused";
  std::vector<int> lsave(5, 0), isave(45, 0);
  std::vector<double> dsave(30, 0.0);
  if (id == 3) {
    write_cell(iw);
    write_cell(iw0);
  }
  double tolerg = 1.0;
  bool geo_ref = (keywrd.find(" GEO_REF") != std::string::npos &&
                  keywrd.find("LOCATE-TS") == std::string::npos);
  bool times = (keywrd.find(" TIMES") != std::string::npos);
  int maxcyc = 100000;
  if (keywrd.find(" CYCLES") != std::string::npos)
    maxcyc = (int)std::lround(reada(keywrd, (int)keywrd.find(" CYCLES")));
  if (keywrd.find("GNORM=") != std::string::npos) {
    tolerg = reada(keywrd, (int)keywrd.find("GNORM="));
    if (keywrd.find(" LET") == std::string::npos && tolerg < 1.e-2) {
      std::printf("  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01\n");
      tolerg = 1.e-2;
    }
  } else {
    tolerg = 1.0;
    if (id != 0) tolerg = id * 2.0 - 1.0;
    if (keywrd.find(" PREC") != std::string::npos) tolerg = tolerg * 0.2;
  }
  double oldstp[13];
  for (int i = 1; i <= 10; ++i) oldstp[i] = 1.0;
  tlast = tleft;
  double tx2 = second(2);
  double tx1 = tx2;
  tleft = tleft - tx2 + time0;
  for (int i = 1; i <= nvar; ++i) nbd[i - 1] = 0;
  bool restrt = (keywrd.find(" RESTART") != std::string::npos);
  int icyc = 0;
  if (restrt) {
    std::fill(isave.begin(), isave.end(), 0);
    std::fill(dsave.begin(), dsave.end(), 0.0);
    nstep = 0;
    double tt0 = 0.0;
    lbfsav(tt0, 0, wa, nwa, iwa, niwa, task, csave, lsave, isave, dsave, nstep, escf);
    time0 = time0 - tt0;
    if (keywrd.find(" 1SCF") != std::string::npos) {
      int i = (int)keywrd.find(" GRAD") + (int)keywrd.find("DERIV");
      compfg(xparam, true, escf, true, grad, i >= 0);
      iflepo = 1;
      for (int k = 1; k <= nvar; ++k) xparam_arr[k] = xparam[k];
      return;
    }
    if (moperr) { for (int k = 1; k <= nvar; ++k) xparam_arr[k] = xparam[k]; return; }
    dcopy1(nvar, &grad[1], 1, &gold[1], 1);
    icyc = nstep;
  } else {
    task = "START";
    nstep = 0;
    icyc = 0;
  }
  double cycmx = 0.0;
  tlast = tleft;
  resfil = false;
  int itry1 = 0;
  double absmin = 1.e6;
  int jcyc = 0;
  double sum = 0.0, const_d = 0.0, tstep = 0.0, rms = 0.0;
  double tprt = 0.0, tt0 = 0.0;
  char txt_ch = ' ';
  double stepmx = 0.0;
  while (true) {
    if (tleft < 1.5 * cycmx || nstep - icyc + 1 > maxcyc) {
      if (nstep - icyc + 1 > maxcyc) {
        std::printf("                    NUMBER OF CYCLES EXCEEDED.  NOW GOING TO FINAL\n");
      } else {
        std::printf("                    THERE IS NOT ENOUGH TIME FOR ANOTHER CYCLE\n                              NOW GOING TO FINAL\n");
      }
      std::printf("\n          - THE CALCULATION IS BEING DUMPED TO DISK\n");
      std::printf("          - RESTART IT USING THE KEYWORD \"RESTART\"\n");
      tt0 = second(1) - time0;
      lbfsav(tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, dsave, nstep, escf);
      iflepo = -1;
      goto l99;
    }
    if (times) timer(" Before SETULB");
    dcopy1(nvar, &xparam[1], 1, &xold[1], 1);
    setulb(nvar, m, &xparam[1], &bot[1], &top[1], &nbd[1], escf, &grad[1], 0.0, 0.0,
           wa, iwa, task, -1, csave, lsave, isave, dsave);
    if (moperr) goto l99;
    if (times) timer(" AFTER SETULB");
    dsave[2] = dsave[2] + 1.e4;
    if (task.substr(0, 2) == "FG") {
      if (nstep > 1) {
        sum = 0.0;
        for (int i = 1; i <= nvar; ++i)
          sum += (xparam[i] - xold[i]) * (xparam[i] - xold[i]);
        sum = std::sqrt(sum);
        int i = std::min(nvar, 10);
        stepmx = std::min(1.0, std::sqrt(ddot1(i, &oldstp[1], 1, &oldstp[1], 1) / i));
        if (sum > stepmx * 2.0) {
          const_d = 2.0 * stepmx / sum;
          for (int i = 1; i <= nvar; ++i)
            xparam[i] = const_d * xparam[i] + (1.0 - const_d) * xold[i];
          sum = 2.0 * stepmx;
        }
        oldstp[nstep % 10 + 1] = sum;
      }
      for (int i = 1; i <= nvar; ++i) {
        if (std::fabs(xparam[i] - xold[i]) > 0.2)
          xparam[i] = xold[i] + std::max(-0.2, std::min(0.2, xparam[i] - xold[i]));
      }
      compfg(xparam, true, escf, true, grad, true);
      if (moperr) goto l99;
      if (absmin - escf < 1.e-7) {
        ++itry1;
        if (itry1 > 900 || (gnorm < 1.0 && itry1 > 9)) {
          std::printf(" HEAT OF FORMATION IS ESSENTIALLY STATIONARY\n");
          iflepo = 3;
          if (best_funct < escf) {
            escf = best_funct;
            gnorm = best_gnorm;
            for (int i = 1; i <= nvar; ++i) xparam[i] = best_xparam[i];
            for (int i = 1; i <= nvar; ++i) grad[i] = best_gradients[i];
            for (int i = 1; i <= natoms; ++i) nc[i] = best_nc[i];
          }
          break;
        }
      } else {
        itry1 = 0;
        absmin = escf;
      }
      if (times) timer(" AFTER COMPFG");
      ++nstep;
      tx2 = second(2);
      tstep = tx2 - tx1;
      cycmx = std::max(tstep, cycmx);
      tx1 = tx2;
      tleft = tleft - tstep;
      if (tlast - tleft > tdump) {
        tlast = tleft;
        resfil = true;
        tt0 = second(1) - time0;
        lbfsav(tt0, 1, wa, nwa, iwa, niwa, task, csave, lsave, isave, dsave, nstep, escf);
        if (moperr) goto l99;
      }
      tleft = std::max(0.0, tleft);
      prttim(tleft, tprt, txt_ch);
      gnorm = std::sqrt(ddot1(nvar, &grad[1], 1, &grad[1], 1));
      if (best_funct > escf) {
        best_gnorm = gnorm;
        best_funct = escf;
        for (int i = 1; i <= nvar; ++i) best_xparam[i] = xparam[i];
        for (int i = 1; i <= nvar; ++i) best_gradients[i] = grad[i];
        for (int i = 1; i <= natoms; ++i) best_nc[i] = nc[i];
      }
      if (id == 3) {
        write_cell(iw);
        write_cell(iw0);
      }
      char txts[2] = {txt_ch, 0};
      char buf[300];
      if (resfil) {
        std::snprintf(buf, sizeof(buf),
                      " RESTART FILE WRITTEN,      TIME LEFT:%6.2f%s  GRAD.:%10.3f HEAT:%14.7g",
                      tprt, txts, std::min(gnorm, 999999.999), escf);
        line = buf;
        std::printf("%s\n", line.c_str());
        to_screen(line);
        if (chanel_C::log) std::printf("%s\n", line.c_str());
        resfil = false;
      } else {
        std::snprintf(buf, sizeof(buf),
                      " CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%s  GRAD.:%10.3f HEAT:%14.7g",
                      nstep, std::min(tstep, 9999.99), tprt, txts,
                      std::min(gnorm, 999999.999), escf);
        if (geo_ref) {
          geo_diff(sum, rms, false);
          char buf2[300];
          std::snprintf(buf2, sizeof(buf2),
                        "   Difference to Geo-Ref:%8.2f = total,%8.4f = Average,%8.4f = RMS movement, in Angstroms",
                        sum, sum / numat, std::sqrt(rms / numat));
          std::printf("%s\n", buf2);
        }
        line = buf;
        std::printf("%s\n", line.c_str());
        if (chanel_C::log) std::printf("%s\n", line.c_str());
        to_screen(line);
      }
      if (nflush != 0 && nstep % nflush == 0) {
        // endfile/backspace on iw/ilog: no-ops
      }
      to_screen("To_file: Geometry optimizing");
      dcopy1(nvar, &grad[1], 1, &gold[1], 1);
      if (gnorm < tolerg) {
        iflepo = 3;
        break;
      }
    } else if (task.substr(0, 5) != "NEW_X") {
      std::printf(" L-BFGS Message:%s\n", task.c_str());
      iflepo = 9;
      break;
    }
  }
  if (gnorm < tolerg)
    std::printf("\n     GRADIENT =%9.5f  IS LESS THAN CUTOFF =%9.5f\n", gnorm, tolerg);
l99:
  if (best_funct + 5.e-4 < escf && keywrd.find(" LET") == std::string::npos) {
    m = 1;
    escf = best_funct;
    gnorm = best_gnorm;
    for (int i = 1; i <= nvar; ++i) xparam[i] = best_xparam[i];
    for (int i = 1; i <= nvar; ++i) grad[i] = best_gradients[i];
    for (int i = 1; i <= natoms; ++i) nc[i] = best_nc[i];
    if (iflepo == -1) {
      for (int i = 1; i <= nvar; ++i) geo[loc[2][i]][loc[1][i]] = xparam[i];
      symtry();
    } else {
      line = " Starting final SCF to re-set LMO's, geometry, etc.";
      std::printf("%s\n", line.c_str());
      to_screen(line);
      last = 1;
      compfg(xparam, true, escf, true, grad, false);
      std::printf("\n          CURRENT BEST VALUE OF HEAT OF FORMATION =%14.6f\n", escf);
      geout(iw);
    }
  } else {
    m = 0;
    if (iflepo != -1 && !mozyme) {
      last = 1;
      compfg(xparam, true, escf, true, grad, false);
    }
  }
  if (iflepo == -1) {
    if (m == 1)
      std::printf("\n          CURRENT BEST VALUE OF HEAT OF FORMATION =%14.6f\n", escf);
    else
      std::printf("\n          CURRENT VALUE OF HEAT OF FORMATION =%14.6f\n", escf);
    if (prt_gradients && keywrd.find(" GRADI") != std::string::npos && mozyme) {
      std::printf("\n\n\n       CURRENT  POINT  AND  DERIVATIVES\n\n");
      prtgra();
    }
    geout(iw);
  }
  for (int k = 1; k <= nvar; ++k) xparam_arr[k] = xparam[k];
  return;
}


// ============================================================================
// lbfsav — save/restore the L-BFGS-B state to a restart file.
// ============================================================================
void lbfsav(double tt0, int mode, std::vector<double>& wa, int nwa,
                   std::vector<int>& iwa, int niwa, std::string& task,
                   std::string& csave, std::vector<int>& lsave, std::vector<int>& isave,
                   std::vector<double>& dsave, int& nstep, double& escf) {
  if (mode == 1) {
    den_in_out(1);
    FILE* fp = std::fopen(restart_fn.c_str(), "wb");
    if (!fp) { mopend("Cannot open restart file for write"); return; }
    std::fwrite(&numat, sizeof(int), 1, fp);
    std::fwrite(&norbs, sizeof(int), 1, fp);
    for (int i = 1; i <= nvar; ++i) { double v = xparam[i]; std::fwrite(&v, sizeof(double), 1, fp); }
    for (int i = 1; i <= nvar; ++i) { double v = grad[i]; std::fwrite(&v, sizeof(double), 1, fp); }
    for (int i = 1; i <= nwa; ++i) { double v = wa[i]; std::fwrite(&v, sizeof(double), 1, fp); }
    for (int i = 1; i <= niwa; ++i) { int v = iwa[i]; std::fwrite(&v, sizeof(int), 1, fp); }
    char tb[61] = {0}, cb[61] = {0};
    std::strncpy(tb, task.c_str(), 60);
    std::strncpy(cb, csave.c_str(), 60);
    std::fwrite(tb, 60, 1, fp);
    std::fwrite(cb, 60, 1, fp);
    for (int i = 1; i <= 4; ++i) { int v = lsave[i]; std::fwrite(&v, sizeof(int), 1, fp); }
    for (int i = 1; i <= 44; ++i) { int v = isave[i]; std::fwrite(&v, sizeof(int), 1, fp); }
    for (int i = 1; i <= 29; ++i) { double v = dsave[i]; std::fwrite(&v, sizeof(double), 1, fp); }
    std::fwrite(&nstep, sizeof(int), 1, fp);
    std::fwrite(&escf, sizeof(double), 1, fp);
    std::fwrite(&nscf, sizeof(int), 1, fp);
    std::fwrite(&tt0, sizeof(double), 1, fp);
    std::fclose(fp);
  } else {
    std::printf("\n          RESTORING DATA FROM DISK\n\n");
    int old_numat = 0, old_norbs = 0;
    FILE* fp = std::fopen(restart_fn.c_str(), "rb");
    int j = 0;
    if (!fp) j = -1;
    if (j == 0 && std::fread(&old_numat, sizeof(int), 1, fp) != 1) j = -1;
    if (j == 0 && std::fread(&old_norbs, sizeof(int), 1, fp) != 1) j = -1;
    for (int i = 1; i <= nvar && j == 0; ++i)
      if (std::fread(&xparam[i], sizeof(double), 1, fp) != 1) j = -1;
    if (j == 0 && std::fread(&grad[1], sizeof(double), (size_t)nvar, fp) != (size_t)nvar) j = -1;
    if (j == 0 && std::fread(&wa[1], sizeof(double), (size_t)nwa, fp) != (size_t)nwa) j = -1;
    if (j == 0 && std::fread(&iwa[1], sizeof(int), (size_t)niwa, fp) != (size_t)niwa) j = -1;
    if (j == 0) {
      char tb[61] = {0}, cb[61] = {0};
      if (std::fread(tb, 60, 1, fp) != 1) j = -1;
      if (std::fread(cb, 60, 1, fp) != 1) j = -1;
      task = tb;
      csave = cb;
      for (int i = 1; i <= 4 && j == 0; ++i)
        if (std::fread(&lsave[i], sizeof(int), 1, fp) != 1) j = -1;
      for (int i = 1; i <= 44 && j == 0; ++i)
        if (std::fread(&isave[i], sizeof(int), 1, fp) != 1) j = -1;
      for (int i = 1; i <= 29 && j == 0; ++i)
        if (std::fread(&dsave[i], sizeof(double), 1, fp) != 1) j = -1;
      if (std::fread(&nstep, sizeof(int), 1, fp) != 1) j = -1;
      if (std::fread(&escf, sizeof(double), 1, fp) != 1) j = -1;
      if (std::fread(&nscf, sizeof(int), 1, fp) != 1) j = -1;
      if (std::fread(&tt0, sizeof(double), 1, fp) != 1) j = -1;
    }
    if (fp) std::fclose(fp);
    if (norbs != old_norbs || numat != old_numat || j != 0) {
      mopend("Restart file read in does not match current data set");
      return;
    }
    int ii = (int)(tt0 / 10000000.0);
    tt0 = tt0 - ii * 10000000;
    std::printf("          TOTAL TIME USED SO FAR:%13.2f SECONDS\n\n", tt0);
  }
}

// ============================================================================
// setulb — L-BFGS-B entry: partition wa and dispatch to mainlb.
// ============================================================================
void setulb(int n, int m, double* x, double* l, double* u, int* nbd,
                   double& f, double* g, double factr, double pgtol,
                   std::vector<double>& wa, std::vector<int>& iwa, std::string& task,
                   int iprint, std::string& csave, std::vector<int>& lsave,
                   std::vector<int>& isave, std::vector<double>& dsave) {
  int ld, lr, lsnd, lss, lsy, lt, lwa, lwn, lws, lwt, lwy, lz;
  if (task == "START") {
    isave[1] = m * n;
    isave[2] = m * m;
    isave[3] = 4 * m * m;
    isave[4] = 1;
    isave[5] = isave[4] + isave[1];
    isave[6] = isave[5] + isave[1];
    isave[7] = isave[6] + isave[2];
    isave[8] = isave[7] + isave[2];
    isave[9] = isave[8];
    isave[10] = isave[9] + isave[2];
    isave[11] = isave[10] + isave[3];
    isave[12] = isave[11] + isave[3];
    isave[13] = isave[12] + n;
    isave[14] = isave[13] + n;
    isave[15] = isave[14] + n;
    isave[16] = isave[15] + n;
  }
  lws = isave[4];
  lwy = isave[5];
  lsy = isave[6];
  lss = isave[7];
  lwt = isave[9];
  lwn = isave[10];
  lsnd = isave[11];
  lz = isave[12];
  lr = isave[13];
  ld = isave[14];
  lt = isave[15];
  lwa = isave[16];
  // Column-major Mat2D views over the wa blocks (1-based offsets).
  Mat2D ws(&wa[lws], n, m);
  Mat2D wy(&wa[lwy], n, m);
  Mat2D sy(&wa[lsy], m, m);
  Mat2D ss(&wa[lss], m, m);
  Mat2D wt(&wa[lwt], m, m);
  Mat2D wn(&wa[lwn], 2 * m, 2 * m);
  Mat2D snd(&wa[lsnd], 2 * m, 2 * m);
  mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, sy, ss, wt, wn, snd,
         &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa], &iwa[1], &iwa[n + 1],
         &iwa[2 * n + 1], task, iprint, csave, lsave, &isave[21], dsave.data());
}


// ============================================================================
// mainlb — L-BFGS-B main kernel (state machine over task).
// ============================================================================
static void mainlb(int n, int m, double* x, double* l, double* u, int* nbd,
                   double& f, double* g, double factr, double pgtol,
                   Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& ss, Mat2D& wt,
                   Mat2D& wn, Mat2D& snd, double* z, double* r, double* d,
                   double* t, double* wa, int* Index, int* iwhere, int* indx2,
                   std::string& task, int iprint, std::string& csave,
                   std::vector<int>& lsave, int* isave, double* dsave) {
  const double one = 1.0, zero = 0.0;
  std::string word = "---";
  bool boxed = false, cnstnd = false, prjctd = false, updatd = false, wrk = false;
  int col = 0, head = 1, i = 0, iback = 0, ifun = 0, ileave = 0, info = 0;
  int itail = 0, iter = 0, itfile = 0, iupdat = 0, iword = 0, k = 0;
  int nact = 0, nenter = 0, nfgv = 0, nfree = 0, nint = 0, nintol = 0, nskip = 0;
  double cachyt = 0, cpu1 = 0, cpu2 = 0, ddum = 0, dnorm = 0, dr = 0, dtd = 0;
  double epsmch = 0, fold = 0, gd = 0, gdold = 0, lnscht = 0, rr = 0;
  double sbgnrm = 0, sbtime = 0, stp = 0, stpmx = 0, theta = 0, time = 0;
  double time1 = 0, time2 = 0, tol = 0, xstep = 0;
  std::vector<int> iorder((size_t)n + 1, 0);
  std::vector<double> wv((size_t)2 * m + 1, 0.0);
  if (task == "START") {
    epsmch = dpmeps();
    fold = 0.0; dnorm = 0.0; time1 = 0.0; time2 = 0.0; cpu1 = 0.0; cpu2 = 0.0;
    gd = 0.0; sbgnrm = 0.0; stp = 0.0; stpmx = 0.0; gdold = 0.0; dtd = 0.0;
    col = 0; head = 1; theta = one; iupdat = 0; updatd = false;
    iback = 0; itail = 0; ifun = 0; iword = 0; nact = 0; ileave = 0; nenter = 0;
    iter = 0; nfgv = 0; nint = 0; nintol = 0; nskip = 0; nfree = n;
    tol = factr * epsmch;
    cachyt = 0; sbtime = 0; lnscht = 0;
    word = "---";
    info = 0;
    itfile = 0;
    if (iprint >= 1) {
      itfile = lbfgs_it;
    }
    errclb(n, m, factr, l, u, nbd, task, info, k);
    if (task.substr(0, 5) == "ERROR") {
      xstep = 0.0;
      k = 0;
      prn3lb(n, x, f, task, iprint, info, itfile, iter, nfgv, nintol, nskip,
             nact, sbgnrm, zero, nint, word, iback, stp, xstep, k, cachyt,
             sbtime, lnscht);
      return;
    } else {
      prn1lb(n, m, l, u, x, iprint, itfile, epsmch);
      active(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, boxed);
      task = "FG_START";
      goto L1300;
    }
  } else {
    prjctd = lsave[1] != 0;
    cnstnd = lsave[2] != 0;
    boxed = lsave[3] != 0;
    updatd = lsave[4] != 0;
    nintol = isave[1];
    itfile = isave[3];
    iback = isave[4];
    nskip = isave[5];
    head = isave[6];
    col = isave[7];
    itail = isave[8];
    iter = isave[9];
    iupdat = isave[10];
    nint = isave[12];
    nfgv = isave[13];
    info = isave[14];
    ifun = isave[15];
    iword = isave[16];
    nfree = isave[17];
    nact = isave[18];
    ileave = isave[19];
    nenter = isave[20];
    theta = dsave[1];
    fold = dsave[2];
    tol = dsave[3];
    dnorm = dsave[4];
    epsmch = dsave[5];
    cpu1 = dsave[6];
    cachyt = dsave[7];
    sbtime = dsave[8];
    lnscht = dsave[9];
    time1 = dsave[10];
    gd = dsave[11];
    stpmx = dsave[12];
    sbgnrm = dsave[13];
    stp = dsave[14];
    gdold = dsave[15];
    dtd = dsave[16];
    cpu2 = 0.0;
    if (task.substr(0, 5) != "FG_LN") {
      if (task.substr(0, 5) == "NEW_X") {
        if (sbgnrm <= pgtol) {
          task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
          goto L1100;
        } else {
          ddum = std::max(std::max(std::fabs(fold), std::fabs(f)), one);
          if ((fold - f) <= tol * ddum) {
            task = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH";
            if (iback >= 10) info = -5;
            goto L1100;
          } else {
            for (i = 1; i <= n; ++i) r[i - 1] = g[i - 1] - r[i - 1];
            rr = ddot1(n, r, 1, r, 1);
            if (std::fabs(stp - one) < 1.e-20) {
              dr = gd - gdold;
              ddum = -gdold;
            } else {
              dr = (gd - gdold) * stp;
              dscal1(n, stp, d, 1);
              ddum = -gdold * stp;
            }
            if (dr <= epsmch * ddum) {
              ++nskip;
              updatd = false;
              if (iprint >= 1)
                std::printf("  ys=%10.3e  -gs=%10.3e BFGS update SKIPPED\n", dr, ddum);
            } else {
              updatd = true;
              ++iupdat;
              matupd(n, m, ws, wy, sy, ss, d, r, itail, iupdat, col, head,
                     theta, rr, dr, stp, dtd);
              formt(m, wt, sy, ss, col, theta, info);
              if (info != 0) {
                if (iprint >= 1)
                  std::printf("\n Nonpositive definiteness in Cholesky factorization in formt;\n   refresh the lbfgs memory and restart the iteration.\n");
                info = 0; col = 0; head = 1; theta = one;
                iupdat = 0; updatd = false;
              }
            }
          }
        }
      } else if (task.substr(0, 5) == "FG_ST") {
        nfgv = 1;
        projgr(n, l, u, nbd, x, g, sbgnrm);
        if (iprint >= 1) {
          std::printf("\nAt iterate%5d    f=%12.5e    |proj g|=%12.5e\n", iter, f, sbgnrm);
          if (itfile > 0)
            std::printf("%4d%4d      -      -      -      -      -%8d%10.3e%10.3e\n",
                        iter, nfgv, 0, sbgnrm, f);
        }
        if (sbgnrm <= pgtol) {
          task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
          goto L1100;
        }
      } else if (task.substr(0, 4) == "STOP") {
        if (task.size() > 6 && task.substr(6, 3) == "CPU") {
          dcopy1(n, t, 1, x, 1);
          dcopy1(n, r, 1, g, 1);
          f = fold;
        }
        goto L1100;
      } else {
        goto L1000;
      }
    // ---- GCP / subspace minimization loop ----
      while (true) {
        if (iprint >= 99) std::printf("\n\nITERATION %5d\n", iter + 1);
        iword = -1;
        if (!cnstnd && col > 0) {
          dcopy1(n, x, 1, z, 1);
          wrk = updatd;
          nint = 0;
        } else {
          cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, z, m, wy, ws, sy, wt,
                 theta, col, head, &wa[1], &wa[2 * m + 1], &wa[4 * m + 1],
                 &wa[6 * m + 1], nint, iprint, sbgnrm, info, epsmch);
          if (info != 0) {
            if (iprint >= 1)
              std::printf("\n Singular triangular system detected;\n   refresh the lbfgs memory and restart the iteration.\n");
            info = 0; col = 0; head = 1; theta = one;
            iupdat = 0; updatd = false;
            cachyt = cachyt + cpu2 - cpu1;
            continue;
          } else {
            cachyt = cachyt + cpu2 - cpu1;
            nintol = nintol + nint;
            freev(n, nfree, Index, nenter, ileave, indx2, iwhere, wrk, updatd,
                  cnstnd, iprint, iter);
            nact = n - nfree;
          }
        }
        if (nfree == 0 || col == 0) break;
        if (wrk)
          formk(n, nfree, Index, nenter, ileave, indx2, iupdat, updatd, wn, snd,
                m, ws, wy, sy, theta, col, head, info);
        if (info != 0) {
          if (iprint >= 1)
            std::printf("\n Nonpositive definiteness in Cholesky factorization in formk;\n   refresh the lbfgs memory and restart the iteration.\n");
          info = 0; col = 0; head = 1; theta = one;
          iupdat = 0; updatd = false;
          sbtime = sbtime + cpu2 - cpu1;
        } else {
          cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, Index, theta, col, head,
                 nfree, cnstnd, info);
          if (info == 0)
            subsm(n, m, nfree, Index, l, u, nbd, z, r, ws, wy, theta, col, head,
                  iword, wv, wn, iprint, info);
          if (info != 0) {
            if (iprint >= 1)
              std::printf("\n Singular triangular system detected;\n   refresh the lbfgs memory and restart the iteration.\n");
            info = 0; col = 0; head = 1; theta = one;
            iupdat = 0; updatd = false;
            sbtime = sbtime + cpu2 - cpu1;
          } else {
            sbtime = sbtime + cpu2 - cpu1;
            break;
          }
        }
      }
      for (i = 1; i <= n; ++i) d[i - 1] = z[i - 1] - x[i - 1];
    }
    // ---- line search loop ----
    while (true) {
      lnsrlb(n, l, u, nbd, x, f, fold, gd, gdold, g, d, r, t, z, stp, dnorm,
             dtd, xstep, stpmx, iter, ifun, iback, nfgv, info, task, boxed,
             cnstnd, csave, isave + 21, dsave + 16);
      if (info == 0 && iback < 20) goto L1200;
      dcopy1(n, t, 1, x, 1);
      dcopy1(n, r, 1, g, 1);
      f = fold;
      if (col == 0) break;
      if (iprint >= 1)
        std::printf("\n Bad direction in the line search;\n   refresh the lbfgs memory and restart the iteration.\n");
      if (info == 0) nfgv = nfgv - 1;
      info = 0; col = 0; head = 1; theta = one;
      iupdat = 0; updatd = false;
      task = "RESTART_FROM_LNSRCH";
      lnscht = lnscht + cpu2 - cpu1;
      while (true) {
        if (iprint >= 99) std::printf("\n\nITERATION %5d\n", iter + 1);
        iword = -1;
        if (!cnstnd && col > 0) {
          dcopy1(n, x, 1, z, 1);
          wrk = updatd;
          nint = 0;
        } else {
          cauchy(n, x, l, u, nbd, g, iorder, iwhere, t, d, z, m, wy, ws, sy, wt,
                 theta, col, head, &wa[1], &wa[2 * m + 1], &wa[4 * m + 1],
                 &wa[6 * m + 1], nint, iprint, sbgnrm, info, epsmch);
          if (info != 0) {
            if (iprint >= 1)
              std::printf("\n Singular triangular system detected;\n   refresh the lbfgs memory and restart the iteration.\n");
            info = 0; col = 0; head = 1; theta = one;
            iupdat = 0; updatd = false;
            cachyt = cachyt + cpu2 - cpu1;
            continue;
          } else {
            cachyt = cachyt + cpu2 - cpu1;
            nintol = nintol + nint;
            freev(n, nfree, Index, nenter, ileave, indx2, iwhere, wrk, updatd,
                  cnstnd, iprint, iter);
            nact = n - nfree;
          }
        }
        if (nfree == 0 || col == 0) break;
        if (wrk)
          formk(n, nfree, Index, nenter, ileave, indx2, iupdat, updatd, wn, snd,
                m, ws, wy, sy, theta, col, head, info);
        if (info != 0) {
          if (iprint >= 1)
            std::printf("\n Nonpositive definiteness in Cholesky factorization in formk;\n   refresh the lbfgs memory and restart the iteration.\n");
          info = 0; col = 0; head = 1; theta = one;
          iupdat = 0; updatd = false;
          sbtime = sbtime + cpu2 - cpu1;
        } else {
          cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, Index, theta, col, head,
                 nfree, cnstnd, info);
          if (info == 0)
            subsm(n, m, nfree, Index, l, u, nbd, z, r, ws, wy, theta, col, head,
                  iword, wv, wn, iprint, info);
          if (info != 0) {
            if (iprint >= 1)
              std::printf("\n Singular triangular system detected;\n   refresh the lbfgs memory and restart the iteration.\n");
            info = 0; col = 0; head = 1; theta = one;
            iupdat = 0; updatd = false;
            sbtime = sbtime + cpu2 - cpu1;
          } else {
            sbtime = sbtime + cpu2 - cpu1;
            break;
          }
        }
      }
      for (i = 1; i <= n; ++i) d[i - 1] = z[i - 1] - x[i - 1];
    }
    if (info == 0) {
      info = -9;
      nfgv = nfgv - 1;
      ifun = ifun - 1;
      iback = iback - 1;
    }
    task = "ABNORMAL_TERMINATION_IN_LNSRCH";
    iter = iter + 1;
    goto L1100;
  L1200:
    if (task.substr(0, 5) != "FG_LN") {
      lnscht = lnscht + cpu2 - cpu1;
      iter = iter + 1;
      projgr(n, l, u, nbd, x, g, sbgnrm);
      prn2lb(n, x, f, g, iprint, itfile, iter, nfgv, nact, sbgnrm, nint, word,
             iword, iback, stp, xstep);
    }
    goto L1300;
  L1100:
    time2 = 0;
    time = time2 - time1;
    k = 0;
    prn3lb(n, x, f, task, iprint, info, itfile, iter, nfgv, nintol, nskip, nact,
           sbgnrm, time, nint, word, iback, stp, xstep, k, cachyt, sbtime, lnscht);
    goto L1300;
  L1000:
    task = "FG_START";
  L1300:
    lsave[1] = prjctd ? 1 : 0;
    lsave[2] = cnstnd ? 1 : 0;
    lsave[3] = boxed ? 1 : 0;
    lsave[4] = updatd ? 1 : 0;
    isave[1] = nintol;
    isave[3] = itfile;
    isave[4] = iback;
    isave[5] = nskip;
    isave[6] = head;
    isave[7] = col;
    isave[8] = itail;
    isave[9] = iter;
    isave[10] = iupdat;
    isave[12] = nint;
    isave[13] = nfgv;
    isave[14] = info;
    isave[15] = ifun;
    isave[16] = iword;
    isave[17] = nfree;
    isave[18] = nact;
    isave[19] = ileave;
    isave[20] = nenter;
    dsave[1] = theta;
    dsave[2] = fold;
    dsave[3] = tol;
    dsave[4] = dnorm;
    dsave[5] = epsmch;
    dsave[6] = cpu1;
    dsave[7] = cachyt;
    dsave[8] = sbtime;
    dsave[9] = lnscht;
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx;
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;
  }
}


// ============================================================================
// active — project x onto the feasible set; classify variables.
// ============================================================================
static void active(int n, double* l, double* u, int* nbd, double* x,
                   int* iwhere, int iprint, bool& prjctd, bool& cnstnd,
                   bool& boxed) {
  const double zero = 0.0;
  int nbdd = 0;
  prjctd = false;
  cnstnd = false;
  boxed = true;
  for (int i = 1; i <= n; ++i) {
    if (nbd[i - 1] > 0) {
      if (nbd[i - 1] <= 2 && x[i - 1] <= l[i - 1]) {
        if (x[i - 1] < l[i - 1]) {
          prjctd = true;
          x[i - 1] = l[i - 1];
        }
        ++nbdd;
      } else if (nbd[i - 1] >= 2 && x[i - 1] >= u[i - 1]) {
        if (x[i - 1] > u[i - 1]) {
          prjctd = true;
          x[i - 1] = u[i - 1];
        }
        ++nbdd;
      }
    }
  }
  for (int i = 1; i <= n; ++i) {
    if (nbd[i - 1] != 2) boxed = false;
    if (nbd[i - 1] == 0) {
      iwhere[i - 1] = -1;
    } else {
      cnstnd = true;
      if (nbd[i - 1] == 2 && u[i - 1] - l[i - 1] <= zero) {
        iwhere[i - 1] = 3;
      } else {
        iwhere[i - 1] = 0;
      }
    }
  }
  if (iprint >= 0) {
    if (prjctd)
      std::printf("The initial X is infeasible.  Restart with its projection.\n");
    if (!cnstnd)
      std::printf("This problem is unconstrained.\n");
  }
  if (iprint > 0)
    std::printf("\nAt X0 %9d variables are exactly at the bounds\n", nbdd);
}

// ============================================================================
// bmv — multiply by the compact representation of B.
// ============================================================================
static void bmv(int m, Mat2D& sy, Mat2D& wt, int col, double* v,
                double* p, int& info) {
  if (col == 0) return;
  p[col + 1 - 1] = v[col + 1 - 1];
  for (int i = 2; i <= col; ++i) {
    int i2 = col + i;
    double sum = 0.0;
    for (int k = 1; k <= i - 1; ++k)
      sum = sum + sy.at(i, k) * v[k - 1] / sy.at(k, k);
    p[i2 - 1] = v[i2 - 1] + sum;
  }
  dtrsl(wt, m, col, &p[col], 11, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i)
    p[i - 1] = v[i - 1] / std::sqrt(sy.at(i, i));
  dtrsl(wt, m, col, &p[col], 1, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i)
    p[i - 1] = -p[i - 1] / std::sqrt(sy.at(i, i));
  for (int i = 1; i <= col; ++i) {
    double sum = 0.0;
    for (int k = i + 1; k <= col; ++k)
      sum = sum + sy.at(k, i) * p[col + k - 1] / sy.at(i, i);
    p[i - 1] = p[i - 1] + sum;
  }
}

// ============================================================================
// cauchy — find the generalized Cauchy point.
// ============================================================================
static void cauchy(int n, double* x, double* l, double* u, int* nbd,
                   const double* g, std::vector<int>& iorder, int* iwhere,
                   double* t, double* d, double* xcp, int m, Mat2D& wy, Mat2D& ws,
                   Mat2D& sy, Mat2D& wt, double& theta, int col, int head,
                   double* p, double* c, double* wbp, double* v, int& nint,
                   int iprint, double sbgnrm, int& info, double epsmch) {
  const double one = 1.0, zero = 0.0;
  if (sbgnrm <= zero) {
    if (iprint >= 0)
      std::printf("Subgnorm = 0.  GCP = X.\n");
    dcopy1(n, x, 1, xcp, 1);
    return;
  }
  bool bnded = true;
  int nfree = n + 1;
  int nbreak = 0;
  int ibkmin = 0;
  double bkmin = zero;
  int col2 = 2 * col;
  double f1 = zero;
  double dibp = zero, dibp2 = zero, dt = zero, wmc = zero, wmp = zero, wmw = zero, zibp = zero;
  if (iprint >= 99)
    std::printf("\n---------------- CAUCHY entered-------------------\n");
  for (int i = 1; i <= col2; ++i) p[i - 1] = zero;
  for (int i = 1; i <= n; ++i) {
    double neggi = -g[i - 1];
    double tl = 0.0, tu = 0.0;
    if (iwhere[i - 1] != 3 && iwhere[i - 1] != -1) {
      if (nbd[i - 1] <= 2) tl = x[i - 1] - l[i - 1];
      if (nbd[i - 1] >= 2) tu = u[i - 1] - x[i - 1];
      bool xlower = nbd[i - 1] <= 2 && tl <= zero;
      bool xupper = nbd[i - 1] >= 2 && tu <= zero;
      iwhere[i - 1] = 0;
      if (xlower) {
        if (neggi <= zero) iwhere[i - 1] = 1;
      } else if (xupper) {
        if (neggi >= zero) iwhere[i - 1] = 2;
      } else if (std::fabs(neggi) <= zero) {
        iwhere[i - 1] = -3;
      }
    }
    int pointr = head;
    if (iwhere[i - 1] != 0 && iwhere[i - 1] != -1) {
      d[i - 1] = zero;
    } else {
      d[i - 1] = neggi;
      f1 = f1 - neggi * neggi;
      for (int j = 1; j <= col; ++j) {
        p[j - 1] = p[j - 1] + wy.at(i, pointr) * neggi;
        p[col + j - 1] = p[col + j - 1] + ws.at(i, pointr) * neggi;
        pointr = pointr % m + 1;
      }
      if (nbd[i - 1] <= 2 && nbd[i - 1] != 0 && neggi < zero) {
        ++nbreak;
        iorder[nbreak] = i;
        t[nbreak - 1] = tl / (-neggi);
        if (nbreak == 1 || t[nbreak - 1] < bkmin) {
          bkmin = t[nbreak - 1];
          ibkmin = nbreak;
        }
      } else if (nbd[i - 1] >= 2 && neggi > zero) {
        ++nbreak;
        iorder[nbreak] = i;
        t[nbreak - 1] = tu / neggi;
        if (nbreak == 1 || t[nbreak - 1] < bkmin) {
          bkmin = t[nbreak - 1];
          ibkmin = nbreak;
        }
      } else {
        --nfree;
        iorder[nfree] = i;
        if (std::fabs(neggi) > zero) bnded = false;
      }
    }
  }
  if (std::fabs(theta - one) > 1.e-20) {
    dscal1(col, theta, &p[col], 1);
  }
  dcopy1(n, x, 1, xcp, 1);
  if (nbreak == 0 && nfree == n + 1) {
    if (iprint > 100) {
      std::printf("Cauchy X =  \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", xcp[i - 1]);
      std::printf("\n");
    }
    return;
  }
  for (int j = 1; j <= col2; ++j) c[j - 1] = zero;
  double f2 = -theta * f1;
  double f2_org = f2;
  if (col > 0) {
    bmv(m, sy, wt, col, p, v, info);
    if (info != 0) return;
    f2 = f2 - ddot1(col2, v, 1, p, 1);
  }
  double dtm = -f1 / f2;
  double tsum = zero;
  nint = 1;
  if (iprint >= 99)
    std::printf("There are %d  breakpoints \n", nbreak);
  if (nbreak != 0) {
    int nleft = nbreak;
    int iter = 1;
    double tj = zero;
    while (true) {
      double tj0 = tj;
      int ibp = 0;
      if (iter == 1) {
        tj = bkmin;
        ibp = iorder[ibkmin];
      } else {
        if (iter == 2) {
          if (ibkmin != nbreak) {
            t[ibkmin - 1] = t[nbreak - 1];
            iorder[ibkmin] = iorder[nbreak];
          }
        }
        hpsolb(nleft, t, iorder, iter - 2);
        tj = t[nleft - 1];
        ibp = iorder[nleft];
      }
      double dt = tj - tj0;
      if (dt != zero && iprint >= 100) {
        std::printf("\nPiece    %3d --f1, f2 at start point %11.4e %11.4e\n",
                    nint, f1, f2);
        std::printf("Distance to the next break point =  %11.4e\n", dt);
        std::printf("Distance to the stationary point =  %11.4e\n", dtm);
      }
      if (dtm < dt) goto L1100;
      tsum = tsum + dt;
      --nleft;
      ++iter;
      double dibp = d[ibp - 1];
      d[ibp - 1] = zero;
      double zibp = 0.0;
      if (dibp > zero) {
        zibp = u[ibp - 1] - x[ibp - 1];
        xcp[ibp - 1] = u[ibp - 1];
        iwhere[ibp - 1] = 2;
      } else {
        zibp = l[ibp - 1] - x[ibp - 1];
        xcp[ibp - 1] = l[ibp - 1];
        iwhere[ibp - 1] = 1;
      }
      if (iprint >= 100)
        std::printf("Variable  %d  is fixed.\n", ibp);
      if (nleft == 0 && nbreak == n) break;
      ++nint;
      double dibp2 = dibp * dibp;
      f1 = f1 + dt * f2 + dibp2 - theta * dibp * zibp;
      f2 = f2 - theta * dibp2;
      if (col > 0) {
        daxpy1(col2, dt, p, 1, c, 1);
        int pointr = head;
        for (int j = 1; j <= col; ++j) {
          wbp[j - 1] = wy.at(ibp, pointr);
          wbp[col + j - 1] = theta * ws.at(ibp, pointr);
          pointr = pointr % m + 1;
        }
        bmv(m, sy, wt, col, wbp, v, info);
        if (info != 0) return;
        double wmc = ddot1(col2, c, 1, v, 1);
        double wmp = ddot1(col2, p, 1, v, 1);
        double wmw = ddot1(col2, wbp, 1, v, 1);
        daxpy1(col2, -dibp, wbp, 1, p, 1);
        f1 = f1 + dibp * wmc;
        f2 = f2 + 2.0 * dibp * wmp - dibp2 * wmw;
      }
      f2 = std::max(epsmch * f2_org, f2);
      if (nleft > 0) {
        dtm = -f1 / f2;
      } else {
        goto L1000;
      }
    }
    dtm = dt;
    goto L1200;
  L1000:
    if (bnded) {
      f1 = zero;
      f2 = zero;
      dtm = zero;
    } else {
      dtm = -f1 / f2;
    }
  }
L1100:
  if (iprint >= 99) {
    std::printf("\n");
    std::printf("GCP found in this segment\n");
    std::printf("Piece    %3d --f1, f2 at start point %11.4e %11.4e\n", nint, f1, f2);
    std::printf("Distance to the stationary point =  %11.4e\n", dtm);
  }
  if (dtm <= zero) dtm = zero;
  tsum = tsum + dtm;
  daxpy1(n, tsum, d, 1, xcp, 1);
L1200:
  if (col > 0) {
    daxpy1(col2, dtm, p, 1, c, 1);
  }
  if (iprint > 100) {
    std::printf("Cauchy X =  \n");
    for (int i = 1; i <= n; ++i) std::printf("%11.4e ", xcp[i - 1]);
    std::printf("\n");
  }
  if (iprint >= 99)
    std::printf("\n---------------- exit CAUCHY----------------------\n");
}

// ============================================================================
// cmprlb — compute r = -Z'B(xcp-x) - Z'g.
// ============================================================================
static void cmprlb(int n, int m, const double* x, const double* g,
                   Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& wt, const double* z,
                   double* r, double* wa, const int* Index, double theta, int col,
                   int head, int nfree, bool cnstnd, int& info) {
  if (!cnstnd && col > 0) {
    for (int i = 1; i <= n; ++i) r[i - 1] = -g[i - 1];
  } else {
    for (int i = 1; i <= nfree; ++i) {
      int k = Index[i - 1];
      r[i - 1] = -theta * (z[k - 1] - x[k - 1]) - g[k - 1];
    }
    bmv(m, sy, wt, col, &wa[2 * m + 1], &wa[1], info);
    if (info != 0) {
      info = -8;
      return;
    }
    int pointr = head;
    for (int j = 1; j <= col; ++j) {
      double a1 = wa[j];
      double a2 = theta * wa[col + j];
      for (int i = 1; i <= nfree; ++i) {
        int k = Index[i - 1];
        r[i - 1] = r[i - 1] + wy.at(k, pointr) * a1 + ws.at(k, pointr) * a2;
      }
      pointr = pointr % m + 1;
    }
  }
}

// ============================================================================
// errclb — check input arguments for errors.
// ============================================================================
static void errclb(int n, int m, double factr, double* l, double* u, int* nbd,
                   std::string& task, int& info, int& k) {
  const double zero = 0.0;
  if (n <= 0) task = "ERROR: N .LE. 0";
  if (m <= 0) task = "ERROR: M .LE. 0";
  if (factr < zero) task = "ERROR: FACTR .LT. 0";
  for (int i = 1; i <= n; ++i) {
    if (nbd[i - 1] < 0 || nbd[i - 1] > 3) {
      task = "ERROR: INVALID NBD";
      info = -6;
      k = i;
    }
    if (nbd[i - 1] == 2) {
      if (l[i - 1] > u[i - 1]) {
        task = "ERROR: NO FEASIBLE SOLUTION";
        info = -7;
        k = i;
      }
    }
  }
}

// ============================================================================
// formk — form and factor the indefinite matrix K = [-D -Y'ZZ'Y/theta  L'-R'].
// ============================================================================
static void formk(int n, int nsub, const int* ind, int nenter, int ileave,
                  const int* indx2, int iupdat, bool updatd, Mat2D& wn,
                  Mat2D& wn1, int m, Mat2D& ws, Mat2D& wy, Mat2D& sy, double theta,
                  int col, int head, int& info) {
  const double zero = 0.0;
  if (updatd) {
    if (iupdat > m) {
      for (int jy = 1; jy <= m - 1; ++jy) {
        int js = m + jy;
        dcopy1(m - jy, &wn1.at(jy + 1, jy + 1), 1, &wn1.at(jy, jy), 1);
        dcopy1(m - jy, &wn1.at(js + 1, js + 1), 1, &wn1.at(js, js), 1);
        dcopy1(m - 1, &wn1.at(m + 2, jy + 1), 1, &wn1.at(m + 1, jy), 1);
      }
    }
    int pbegin = 1;
    int pend = nsub;
    int dbegin = nsub + 1;
    int dend = n;
    int iy = col;
    int is = m + col;
    int ipntr = head + col - 1;
    if (ipntr > m) ipntr = ipntr - m;
    int jpntr = head;
    for (int jy = 1; jy <= col; ++jy) {
      int js = m + jy;
      double temp1 = zero, temp2 = zero, temp3 = zero;
      for (int k = pbegin; k <= pend; ++k) {
        int k1 = ind[k - 1];
        temp1 = temp1 + wy.at(k1, ipntr) * wy.at(k1, jpntr);
      }
      for (int k = dbegin; k <= dend; ++k) {
        int k1 = ind[k - 1];
        temp2 = temp2 + ws.at(k1, ipntr) * ws.at(k1, jpntr);
        temp3 = temp3 + ws.at(k1, ipntr) * wy.at(k1, jpntr);
      }
      wn1.at(iy, jy) = temp1;
      wn1.at(is, js) = temp2;
      wn1.at(is, jy) = temp3;
      jpntr = jpntr % m + 1;
    }
    int jy = col;
    jpntr = head + col - 1;
    if (jpntr > m) jpntr = jpntr - m;
    ipntr = head;
    for (int i = 1; i <= col; ++i) {
      int is2 = m + i;
      double temp3 = zero;
      for (int k = pbegin; k <= pend; ++k) {
        int k1 = ind[k - 1];
        temp3 = temp3 + ws.at(k1, ipntr) * wy.at(k1, jpntr);
      }
      ipntr = ipntr % m + 1;
      wn1.at(is2, jy) = temp3;
    }
  }
  int upcl = updatd ? col - 1 : col;
  int ipntr = head;
  for (int iy = 1; iy <= upcl; ++iy) {
    int is = m + iy;
    int jpntr = head;
    for (int jy = 1; jy <= iy; ++jy) {
      int js = m + jy;
      double temp1 = zero, temp2 = zero, temp3 = zero, temp4 = zero;
      for (int k = 1; k <= nenter; ++k) {
        int k1 = indx2[k - 1];
        temp1 = temp1 + wy.at(k1, ipntr) * wy.at(k1, jpntr);
        temp2 = temp2 + ws.at(k1, ipntr) * ws.at(k1, jpntr);
      }
      for (int k = ileave; k <= n; ++k) {
        int k1 = indx2[k - 1];
        temp3 = temp3 + wy.at(k1, ipntr) * wy.at(k1, jpntr);
        temp4 = temp4 + ws.at(k1, ipntr) * ws.at(k1, jpntr);
      }
      wn1.at(iy, jy) = wn1.at(iy, jy) + temp1 - temp3;
      wn1.at(is, js) = wn1.at(is, js) - temp2 + temp4;
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }
  ipntr = head;
  for (int is = m + 1; is <= m + upcl; ++is) {
    int jpntr = head;
    for (int jy = 1; jy <= upcl; ++jy) {
      double temp1 = zero, temp3 = zero;
      for (int k = 1; k <= nenter; ++k) {
        int k1 = indx2[k - 1];
        temp1 = temp1 + ws.at(k1, ipntr) * wy.at(k1, jpntr);
      }
      for (int k = ileave; k <= n; ++k) {
        int k1 = indx2[k - 1];
        temp3 = temp3 + ws.at(k1, ipntr) * wy.at(k1, jpntr);
      }
      if (is <= jy + m)
        wn1.at(is, jy) = wn1.at(is, jy) + temp1 - temp3;
      else
        wn1.at(is, jy) = wn1.at(is, jy) - temp1 + temp3;
      jpntr = jpntr % m + 1;
    }
    ipntr = ipntr % m + 1;
  }
  int m2 = 2 * m;
  for (int iy = 1; iy <= col; ++iy) {
    int is = col + iy;
    int is1 = m + iy;
    for (int jy = 1; jy <= iy; ++jy) {
      int js = col + jy;
      int js1 = m + jy;
      wn.at(jy, iy) = wn1.at(iy, jy) / theta;
      wn.at(js, is) = wn1.at(is1, js1) * theta;
    }
    for (int jy = 1; jy <= iy - 1; ++jy)
      wn.at(jy, is) = -wn1.at(is1, jy);
    for (int jy = iy; jy <= col; ++jy)
      wn.at(jy, is) = wn1.at(is1, jy);
    wn.at(iy, iy) = wn.at(iy, iy) + sy.at(iy, iy);
  }
  dpofa(wn, m2, col, info);
  if (info != 0) {
    info = -1;
    return;
  }
  int col2 = 2 * col;
  for (int js = col + 1; js <= col2; ++js)
    dtrsl(wn, m2, col, &wn.at(1, js), 11, info);
  for (int is = col + 1; is <= col2; ++is)
    for (int js = is; js <= col2; ++js)
      wn.at(is, js) = wn.at(is, js) + ddot1(col, &wn.at(1, is), 1, &wn.at(1, js), 1);
  dpofa(wn, m2, col, info, col + 1, col + 1);
  if (info != 0) info = -2;
}

// ============================================================================
// formt — form T = theta*SS + L*D^(-1)*L' and its Cholesky factor.
// ============================================================================
static void formt(int m, Mat2D& wt, Mat2D& sy, Mat2D& ss, int col, double theta, int& info) {
  const double zero = 0.0;
  for (int j = 1; j <= col; ++j)
    wt.at(1, j) = theta * ss.at(1, j);
  for (int i = 2; i <= col; ++i)
    for (int j = i; j <= col; ++j) {
      int k1 = std::min(i, j) - 1;
      double ddum = zero;
      for (int k = 1; k <= k1; ++k)
        ddum = ddum + sy.at(i, k) * sy.at(j, k) / sy.at(k, k);
      wt.at(i, j) = ddum + theta * ss.at(i, j);
    }
  dpofa(wt, m, col, info);
  if (info != 0) info = -3;
}

// ============================================================================
// freev — count entering/leaving variables; re-index free set.
// ============================================================================
static void freev(int n, int& nfree, int* Index, int& nenter, int& ileave,
                  int* indx2, const int* iwhere, bool& wrk,
                  bool updatd, bool cnstnd, int iprint, int iter) {
  nenter = 0;
  ileave = n + 1;
  if (iter > 0 && cnstnd) {
    for (int i = 1; i <= nfree; ++i) {
      int k = Index[i - 1];
      if (iwhere[k - 1] > 0) {
        --ileave;
        indx2[ileave - 1] = k;
        if (iprint >= 100)
          std::printf("Variable %d leaves the set of free variables\n", k);
      }
    }
    for (int i = 1 + nfree; i <= n; ++i) {
      int k = Index[i - 1];
      if (iwhere[k - 1] <= 0) {
        ++nenter;
        indx2[nenter - 1] = k;
        if (iprint >= 100)
          std::printf("Variable %d enters the set of free variables\n", k);
      }
    }
    if (iprint >= 99)
      std::printf("%d variables leave; %d variables enter\n", n + 1 - ileave, nenter);
  }
  wrk = (ileave < n + 1) || (nenter > 0) || updatd;
  nfree = 0;
  int iact = n + 1;
  for (int i = 1; i <= n; ++i) {
    if (iwhere[i - 1] <= 0) {
      ++nfree;
      Index[nfree - 1] = i;
    } else {
      --iact;
      Index[iact - 1] = i;
    }
  }
  if (iprint >= 99)
    std::printf("%d variables are free at GCP %d\n", nfree, iter + 1);
}

// ============================================================================
// hpsolb — heap sort the breakpoints (used by cauchy).
// ============================================================================
static void hpsolb(int n, double* t, std::vector<int>& iorder, int iheap) {
  if (iheap == 0) {
    for (int k = 2; k <= n; ++k) {
      double ddum = t[k - 1];
      int indxin = iorder[k];
      int i = k;
      while (i > 1) {
        int j = i / 2;
        if (ddum < t[j - 1]) {
          t[i - 1] = t[j - 1];
          iorder[i] = iorder[j];
          i = j;
        } else {
          break;
        }
      }
      t[i - 1] = ddum;
      iorder[i] = indxin;
    }
  }
  if (n <= 1) return;
  int i = 1;
  double out = t[1 - 1];
  int indxou = iorder[1];
  double ddum = t[n - 1];
  int indxin = iorder[n];
  while (true) {
    int j = i + i;
    if (j > n - 1) break;
    if (t[j + 1 - 1] < t[j - 1]) j = j + 1;
    if (t[j - 1] < ddum) {
      t[i - 1] = t[j - 1];
      iorder[i] = iorder[j];
      i = j;
    } else {
      break;
    }
  }
  t[i - 1] = ddum;
  iorder[i] = indxin;
  t[n - 1] = out;
  iorder[n] = indxou;
}

// ============================================================================
// lnsrlb — line search along d (drives dcsrch).
// ============================================================================
static void lnsrlb(int n, double* l, double* u, int* nbd, double* x, double f,
                   double& fold, double& gd, double& gdold, double* g,
                   double* d, double* r, double* t, double* z, double& stp,
                   double& dnorm, double& dtd, double& xstep, double& stpmx,
                   int iter, int& ifun, int& iback, int& nfgv, int& info,
                   std::string& task, bool boxed, bool cnstnd, std::string& csave,
                   int* isave, double* dsave) {
  const double one = 1.0, zero = 0.0;
  const double big = 1.e5;
  const double ftol = 1.0e-3, gtol = 0.9, xtol = 0.1;
  if (task.substr(0, 5) != "FG_LN") {
    dtd = ddot1(n, d, 1, d, 1);
    dnorm = std::sqrt(dtd);
    stpmx = big;
    if (cnstnd) {
      if (iter == 0) {
        stpmx = one;
      } else {
        for (int i = 1; i <= n; ++i) {
          double a1 = d[i - 1];
          if (nbd[i - 1] != 0) {
            if (a1 < zero && nbd[i - 1] <= 2) {
              double a2 = l[i - 1] - x[i - 1];
              if (a2 >= zero) stpmx = zero;
              else if (a1 * stpmx < a2) stpmx = a2 / a1;
            } else if (a1 > zero && nbd[i - 1] >= 2) {
              double a2 = u[i - 1] - x[i - 1];
              if (a2 <= zero) stpmx = zero;
              else if (a1 * stpmx > a2) stpmx = a2 / a1;
            }
          }
        }
      }
    }
    if (iter == 0 && !boxed) stp = std::min(one / dnorm, stpmx);
    else stp = one;
    dcopy1(n, x, 1, t, 1);
    dcopy1(n, g, 1, r, 1);
    fold = f;
    ifun = 0;
    iback = 0;
    csave = "START";
  }
  gd = ddot1(n, g, 1, d, 1);
  if (ifun == 0) {
    gdold = gd;
    if (gd >= zero) {
      info = -4;
      return;
    }
  }
  dcsrch(f, gd, stp, ftol, gtol, xtol, zero, stpmx, csave, isave, dsave);
  xstep = stp * dnorm;
  if (csave.substr(0, 4) == "CONV" || csave.substr(0, 4) == "WARN") {
    task = "NEW_X";
    return;
  }
  task = "FG_LNSRCH";
  ++ifun;
  ++nfgv;
  iback = ifun - 1;
  if (std::fabs(stp - one) < 1.e-20) {
    dcopy1(n, z, 1, x, 1);
  } else {
    for (int i = 1; i <= n; ++i)
      x[i - 1] = stp * d[i - 1] + t[i - 1];
  }
}

// ============================================================================
// matupd — update the L-BFGS matrices WS, WY, SY, SS.
// ============================================================================
static void matupd(int n, int m, Mat2D& ws, Mat2D& wy, Mat2D& sy, Mat2D& ss,
                   double* d, double* r, int& itail,
                   int iupdat, int& col, int& head, double& theta, double rr,
                   double dr, double stp, double dtd) {
  const double one = 1.0;
  if (iupdat <= m) {
    col = iupdat;
    itail = (head + iupdat - 2) % m + 1;
  } else {
    itail = itail % m + 1;
    head = head % m + 1;
  }
  dcopy1(n, d, 1, ws.col(itail), 1);
  dcopy1(n, r, 1, wy.col(itail), 1);
  theta = rr / dr;
  if (iupdat > m) {
    for (int j = 1; j <= col - 1; ++j) {
      dcopy1(j, &ss.at(2, j + 1), 1, &ss.at(1, j), 1);
      dcopy1(col - j, &sy.at(j + 1, j + 1), 1, &sy.at(j, j), 1);
    }
  }
  int pointr = head;
  for (int j = 1; j <= col - 1; ++j) {
    sy.at(col, j) = ddot1(n, d, 1, wy.col(pointr), 1);
    ss.at(j, col) = ddot1(n, ws.col(pointr), 1, d, 1);
    pointr = pointr % m + 1;
  }
  if (std::fabs(stp - one) < 1.e-20)
    ss.at(col, col) = dtd;
  else
    ss.at(col, col) = stp * stp * dtd;
  sy.at(col, col) = dr;
}


// ============================================================================
// prn1lb — print banner / initial info of L-BFGS-B.
// ============================================================================
static void prn1lb(int n, int m, double* l, double* u, double* x, int iprint,
                   int itfile, double epsmch) {
  (void)l; (void)u; (void)x;
  if (iprint < 0) return;
  std::printf("RUNNING THE L-BFGS-B CODE\n\n           * * *\n\nMachine precision =%10.3e\n",
              epsmch);
  std::printf("N = %d    M = %d\n", n, m);
  if (iprint >= 1) {
    if (itfile > 0) {
      std::printf("RUNNING THE L-BFGS-B CODE\n\nit    = iteration number\n");
      std::printf("nf    = number of function evaluations\n");
      std::printf("nint  = number of segments explored during the Cauchy search\n");
      std::printf("nact  = number of active bounds at the generalized Cauchy point\n");
      std::printf("sub   = manner in which the subspace minimization terminated:\n");
      std::printf("        con = converged, bnd = a bound was reached\n");
      std::printf("itls  = number of iterations performed in the line search\n");
      std::printf("stepl = step length used\n");
      std::printf("tstep = norm of the displacement (total step)\n");
      std::printf("projg = norm of the projected gradient\n");
      std::printf("f     = function value\n\n           * * *\n\nMachine precision =%10.3e\n",
                  epsmch);
      std::printf("N = %d    M = %d\n", n, m);
      std::printf("\n   it   nf  nint  nact  sub  itls   stepl    tstep      projg         f\n");
    }
    if (iprint > 100) {
      std::printf("\nL =   \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", l[i - 1]);
      std::printf("\nX0 =  \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", x[i - 1]);
      std::printf("\nU =   \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", u[i - 1]);
      std::printf("\n");
    }
  }
}

// ============================================================================
// prn2lb — print info about a new iterate.
// ============================================================================
static void prn2lb(int n, double* x, double f, double* g, int iprint, int itfile,
                   int iter, int nfgv, int nact, double sbgnrm, int nint,
                   std::string& word, int iword, int iback, double stp, double xstep) {
  (void)nint;
  if (iword == 0) word = "con";
  else if (iword == 1) word = "bnd";
  else if (iword == 5) word = "TNT";
  else word = "---";
  if (iprint >= 99) {
    std::printf("LINE SEARCH %d times; norm of step = %g\n", iback, xstep);
    std::printf("\nAt iterate%5d    f=%12.5e    |proj g|=%12.5e\n", iter, f, sbgnrm);
    if (iprint > 100) {
      std::printf("\nX =   \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", x[i - 1]);
      std::printf("\nG =   \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", g[i - 1]);
      std::printf("\n");
    }
  } else if (iprint > 0) {
    if (iter % iprint == 0)
      std::printf("\nAt iterate%5d    f=%12.5e    |proj g|=%12.5e\n", iter, f, sbgnrm);
  }
  if (iprint >= 1 && itfile > 0)
    std::printf("%4d%4d%5d%5d  %s%4d%8.1e%8.1e%10.3e%10.3e\n",
                iter, nfgv, nint, nact, word.c_str(), iback, stp, xstep, sbgnrm, f);
}

// ============================================================================
// prn3lb — print termination info of L-BFGS-B.
// ============================================================================
static void prn3lb(int n, double* x, double f, const std::string& task, int iprint,
                   int info, int itfile, int iter, int nfgv, int nintol, int nskip,
                   int nact, double sbgnrm, double time, int nint,
                   const std::string& word, int iback, double stp, double xstep, int k,
                   double cachyt, double sbtime, double lnscht) {
  (void)nint; (void)stp; (void)xstep;
  if (task.substr(0, 5) != "ERROR") {
    if (iprint >= 0) {
      std::printf("\n           * * *\n\nTit   = total number of iterations\n");
      std::printf("Tnf   = total number of function evaluations\n");
      std::printf("Tnint = total number of segments explored during Cauchy searches\n");
      std::printf("Skip  = number of BFGS updates skipped\n");
      std::printf("Nact  = number of active bounds at final generalized Cauchy point\n");
      std::printf("Projg = norm of the final projected gradient\n");
      std::printf("F     = final function value\n\n           * * *\n");
      std::printf("\n   N  Tit  Tnf  Tnint  Skip  Nact       Projg           F\n");
      std::printf("%5d%4d%4d%7d%6d%6d%10.3e%10.3e\n",
                  n, iter, nfgv, nintol, nskip, nact, sbgnrm, f);
      if (iprint >= 100) {
        std::printf("\nX =   \n");
        for (int i = 1; i <= n; ++i) std::printf("%11.4e ", x[i - 1]);
        std::printf("\n");
      }
      if (iprint >= 1) std::printf(" F =%g\n", f);
    }
  }
  if (iprint < 0) return;
  std::printf("\n%s\n", task.c_str());
  if (info != 0) {
    if (info == -1)
      std::printf("\n Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n");
    if (info == -2)
      std::printf("\n Matrix in 2st Cholesky factorization in formk is not Pos. Def.\n");
    if (info == -3)
      std::printf("\n Matrix in the Cholesky factorization in formt is not Pos. Def.\n");
    if (info == -4)
      std::printf("\n Derivative >= 0, backtracking line search impossible.\n   Previous x, f and g restored.\n Possible causes: 1 error in function or gradient evaluation;\n                 2 rounding errors dominate computation.\n");
    if (info == -5)
      std::printf("\n Warning:  more than 10 function and gradient\n   evaluations in the last line search.  Termination\n   may possibly be caused by a bad search direction.\n");
    if (info == -6)
      std::printf(" Input nbd(%d) is invalid.\n", k);
    if (info == -7)
      std::printf(" l(%d) > u(%d).  No feasible solution.\n", k, k);
    if (info == -8)
      std::printf("\n The triangular system is singular.\n");
    if (info == -9)
      std::printf("\n Line search cannot locate an adequate point after 20 function\n  and gradient evaluations.  Previous x, f and g restored.\n Possible causes: 1 error in function or gradient evaluation;\n                 2 rounding error dominate computation.\n");
  }
  if (iprint >= 1) {
    std::printf("\n Cauchy                time%10.3e seconds.\n", cachyt);
    std::printf(" Subspace minimization time%10.3e seconds.\n", sbtime);
    std::printf(" Line search           time%10.3e seconds.\n", lnscht);
  }
  std::printf("\n Total User time%10.3e seconds.\n\n", time);
  if (iprint < 1) return;
  if (info == -4 || info == -9) {
    if (itfile > 0)
      std::printf("%4d%4d%5d%5d  %s%4d%8.1e%8.1e      -          -\n",
                  iter, nfgv, nint, nact, word.c_str(), iback, stp, xstep);
  }
  if (itfile > 0) {
    std::printf("\n%s\n", task.c_str());
    if (info == -1) std::printf("\n Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n");
    if (info == -2) std::printf("\n Matrix in 2st Cholesky factorization in formk is not Pos. Def.\n");
    if (info == -3) std::printf("\n Matrix in the Cholesky factorization in formt is not Pos. Def.\n");
    if (info == -4)
      std::printf("\n Derivative >= 0, backtracking line search impossible.\n   Previous x, f and g restored.\n Possible causes: 1 error in function or gradient evaluation;\n                 2 rounding errors dominate computation.\n");
    if (info == -5)
      std::printf("\n Warning:  more than 10 function and gradient\n   evaluations in the last line search.  Termination\n   may possibly be caused by a bad search direction.\n");
    if (info == -8) std::printf("\n The triangular system is singular.\n");
    if (info == -9)
      std::printf("\n Line search cannot locate an adequate point after 20 function\n  and gradient evaluations.  Previous x, f and g restored.\n Possible causes: 1 error in function or gradient evaluation;\n                 2 rounding error dominate computation.\n");
    std::printf("\n Total User time%10.3e seconds.\n\n", time);
  }
}

// ============================================================================
// projgr — infinity norm of the projected gradient.
// ============================================================================
static void projgr(int n, double* l, double* u, int* nbd, double* x, double* g,
                   double& sbgnrm) {
  const double zero = 0.0;
  sbgnrm = zero;
  for (int i = 1; i <= n; ++i) {
    double gi = g[i - 1];
    if (nbd[i - 1] != 0) {
      if (gi < zero) {
        if (nbd[i - 1] >= 2)
          gi = std::max((x[i - 1] - u[i - 1]), gi);
      } else if (nbd[i - 1] <= 2) {
        gi = std::min((x[i - 1] - l[i - 1]), gi);
      }
    }
    sbgnrm = std::max(sbgnrm, std::fabs(gi));
  }
}

// ============================================================================
// subsm — subspace minimization (direct method).
// ============================================================================
static void subsm(int n, int m, int nsub, const int* ind, double* l, double* u,
                  int* nbd, double* x, double* d, Mat2D& ws, Mat2D& wy,
                  double theta, int col, int head, int& iword, std::vector<double>& wv,
                  Mat2D& wn, int iprint, int& info) {
  const double one = 1.0, zero = 0.0;
  (void)one;
  if (nsub <= 0) return;
  if (iprint >= 99)
    std::printf("\n----------------SUBSM entered-----------------\n");
  int pointr = head;
  for (int i = 1; i <= col; ++i) {
    double temp1 = zero, temp2 = zero;
    for (int j = 1; j <= nsub; ++j) {
      int k = ind[j - 1];
      temp1 = temp1 + wy.at(k, pointr) * d[j - 1];
      temp2 = temp2 + ws.at(k, pointr) * d[j - 1];
    }
    wv[i] = temp1;
    wv[col + i] = theta * temp2;
    pointr = pointr % m + 1;
  }
  int m2 = 2 * m;
  int col2 = 2 * col;
  dtrsl(wn, m2, col2, &wv[1], 11, info);
  if (info != 0) return;
  for (int i = 1; i <= col; ++i) wv[i] = -wv[i];
  dtrsl(wn, m2, col2, &wv[1], 1, info);
  if (info != 0) return;
  pointr = head;
  for (int jy = 1; jy <= col; ++jy) {
    int js = col + jy;
    for (int i = 1; i <= nsub; ++i) {
      int k = ind[i - 1];
      d[i - 1] = d[i - 1] + wy.at(k, pointr) * wv[jy] / theta + ws.at(k, pointr) * wv[js];
    }
    pointr = pointr % m + 1;
  }
  for (int i = 1; i <= nsub; ++i) d[i - 1] = d[i - 1] / theta;
  double alpha = one;
  double temp1 = alpha;
  int ibd = 0;
  for (int i = 1; i <= nsub; ++i) {
    int k = ind[i - 1];
    double dk = d[i - 1];
    if (nbd[k - 1] != 0) {
      if (dk < zero && nbd[k - 1] <= 2) {
        double temp2 = l[k - 1] - x[k - 1];
        if (temp2 >= zero) temp1 = zero;
        else if (dk * alpha < temp2) temp1 = temp2 / dk;
      } else if (dk > zero && nbd[k - 1] >= 2) {
        double temp2 = u[k - 1] - x[k - 1];
        if (temp2 <= zero) temp1 = zero;
        else if (dk * alpha > temp2) temp1 = temp2 / dk;
      }
      if (temp1 < alpha) {
        alpha = temp1;
        ibd = i;
      }
    }
  }
  if (alpha < one) {
    double dk = d[ibd - 1];
    int k = ind[ibd - 1];
    if (dk > zero) {
      x[k - 1] = u[k - 1];
      d[ibd - 1] = zero;
    } else if (dk < zero) {
      x[k - 1] = l[k - 1];
      d[ibd - 1] = zero;
    }
  }
  for (int i = 1; i <= nsub; ++i) {
    int k = ind[i - 1];
    x[k - 1] = x[k - 1] + alpha * d[i - 1];
  }
  if (iprint >= 99) {
    if (alpha < one)
      std::printf("ALPHA = %8.5f backtrack to the BOX\n", alpha);
    else
      std::printf("SM solution inside the box\n");
    if (iprint > 100) {
      std::printf("Subspace solution X =  \n");
      for (int i = 1; i <= n; ++i) std::printf("%11.4e ", x[i - 1]);
      std::printf("\n");
    }
  }
  iword = (alpha < one) ? 1 : 0;
  if (iprint >= 99)
    std::printf("\n----------------exit SUBSM --------------------\n");
}


// ============================================================================
// dcsrch — More-Thuente line search (task protocol START/FG).
// ============================================================================
static void dcsrch(double f, double g, double& stp, double ftol, double gtol,
                   double xtol, double stpmin, double stpmax, std::string& task,
                   int* isave, double* dsave) {
  const double zero = 0.0, p5 = 0.5, p66 = 0.66, xtrapl = 1.1, xtrapu = 4.0;
  bool brackt;
  int stage;
  double finit, fm, ftest, fx, fxm, fy, fym, ginit, gm, gtest, gx, gxm, gy, gym;
  double stmax, stmin, stx, sty, width, width1;
  if (task.substr(0, 5) == "START") {
    if (stp < stpmin) task = "ERROR: STP .LT. STPMIN";
    if (stp > stpmax) task = "ERROR: STP .GT. STPMAX";
    if (g >= zero) task = "ERROR: INITIAL G .GE. ZERO";
    if (ftol < zero) task = "ERROR: FTOL .LT. ZERO";
    if (gtol < zero) task = "ERROR: GTOL .LT. ZERO";
    if (xtol < zero) task = "ERROR: XTOL .LT. ZERO";
    if (stpmin < zero) task = "ERROR: STPMIN .LT. ZERO";
    if (stpmax < stpmin) task = "ERROR: STPMAX .LT. STPMIN";
    if (task.substr(0, 5) == "ERROR") return;
    brackt = false;
    stage = 1;
    finit = f;
    ginit = g;
    gtest = ftol * ginit;
    width = stpmax - stpmin;
    width1 = width / p5;
    stx = 1.e-5;
    fx = finit;
    gx = ginit;
    sty = 1.e-4;
    fy = finit;
    gy = ginit;
    stmin = zero;
    stmax = stp + xtrapu * stp;
    task = "FG";
  } else {
    brackt = isave[1] == 1;
    stage = isave[2];
    ginit = dsave[1];
    gtest = dsave[2];
    gx = dsave[3];
    gy = dsave[4];
    finit = dsave[5];
    fx = dsave[6];
    fy = dsave[7];
    stx = dsave[8];
    sty = dsave[9];
    stmin = dsave[10];
    stmax = dsave[11];
    width = dsave[12];
    width1 = dsave[13];
    ftest = finit + stp * gtest;
    if (stage == 1 && f <= ftest && g >= zero) stage = 2;
    if (brackt && (stp <= stmin || stp >= stmax))
      task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS";
    if (brackt && stmax - stmin <= xtol * stmax)
      task = "WARNING: XTOL TEST SATISFIED";
    if (std::fabs(stp - stpmax) < 1.e-20 && f <= ftest && g <= gtest)
      task = "WARNING: STP = STPMAX";
    if (std::fabs(stp - stpmin) < 1.e-20 && (f > ftest || g >= gtest))
      task = "WARNING: STP = STPMIN";
    if (f <= ftest && std::fabs(g) <= gtol * (-ginit))
      task = "CONVERGENCE";
    if (task.substr(0, 4) != "WARN" && task.substr(0, 4) != "CONV") {
      if (stage == 1 && f <= fx && f > ftest) {
        fm = f - stp * gtest;
        fxm = fx - stx * gtest;
        fym = fy - sty * gtest;
        gm = g - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;
        dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
      } else {
        dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
      }
      if (brackt) {
        if (std::fabs(sty - stx) >= p66 * width1) stp = stx + p5 * (sty - stx);
        width1 = width;
        width = std::fabs(sty - stx);
      }
      if (brackt) {
        stmin = std::min(stx, sty);
        stmax = std::max(stx, sty);
      } else {
        stmin = stp + xtrapl * (stp - stx);
        stmax = stp + xtrapu * (stp - stx);
      }
      stp = std::max(stp, stpmin);
      stp = std::min(stp, stpmax);
      if ((brackt && (stp <= stmin || stp >= stmax)) ||
          (brackt && stmax - stmin <= xtol * stmax)) {
        stp = stx;
      }
      task = "FG";
    }
  }
  isave[1] = brackt ? 1 : 0;
  isave[2] = stage;
  dsave[1] = ginit;
  dsave[2] = gtest;
  dsave[3] = gx;
  dsave[4] = gy;
  dsave[5] = finit;
  dsave[6] = fx;
  dsave[7] = fy;
  dsave[8] = stx;
  dsave[9] = sty;
  dsave[10] = stmin;
  dsave[11] = stmax;
  dsave[12] = width;
  dsave[13] = width1;
}

// ============================================================================
// dcstep — one step of the More-Thuente line search.
// ============================================================================
static void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy,
                   double& dy, double& stp, double fp, double dp, bool& brackt,
                   double stpmin, double stpmax) {
  const double zero = 0.0, p66 = 0.66, two = 2.0, three = 3.0;
  double sgnd = dp * (dx / std::fabs(dx));
  if (std::fabs(stp - stx) < 1.e-5) stp = stp + 1.e-5;
  double stpf = 0.0;
  if (fp > fx) {
    double theta = three * (fx - fp) / (stp - stx) + dx + dp;
    double s = std::max(std::max(std::fabs(theta), std::fabs(dx)), std::fabs(dp));
    double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp < stx) gamma = -gamma;
    double p = (gamma - dx) + theta;
    double q = ((gamma - dx) + gamma) + dp;
    double r = p / q;
    double stpc = stx + r * (stp - stx);
    double stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / two) * (stp - stx);
    if (std::fabs(stpc - stx) < std::fabs(stpq - stx))
      stpf = stpc;
    else
      stpf = stpc + (stpq - stpc) / two;
    brackt = true;
  } else if (sgnd < zero) {
    double theta = three * (fx - fp) / (stp - stx) + dx + dp;
    double s = std::max(std::max(std::fabs(theta), std::fabs(dx)), std::fabs(dp));
    double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
    if (stp > stx) gamma = -gamma;
    double p = (gamma - dp) + theta;
    double q = ((gamma - dp) + gamma) + dx;
    double r = p / q;
    double stpc = stp + r * (stx - stp);
    double stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (std::fabs(stpc - stp) > std::fabs(stpq - stp))
      stpf = stpc;
    else
      stpf = stpq;
    brackt = true;
  } else if (std::fabs(dp) < std::fabs(dx)) {
    double theta = three * (fx - fp) / (stp - stx) + dx + dp;
    double s = std::max(std::max(std::fabs(theta), std::fabs(dx)), std::fabs(dp));
    double gamma = s * std::sqrt(std::max(zero, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
    if (stp > stx) gamma = -gamma;
    double p = (gamma - dp) + theta;
    double q = (gamma + (dx - dp)) + gamma;
    double r = p / q;
    double stpc = 0.0;
    if (r < zero && gamma != zero)
      stpc = stp + r * (stx - stp);
    else if (stp > stx)
      stpc = stpmax;
    else
      stpc = stpmin;
    double stpq = stp + (dp / (dp - dx)) * (stx - stp);
    if (brackt) {
      if (std::fabs(stpc - stp) < std::fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      if (stp > stx)
        stpf = std::min(stp + p66 * (sty - stp), stpf);
      else
        stpf = std::max(stp + p66 * (sty - stp), stpf);
    } else {
      if (std::fabs(stpc - stp) > std::fabs(stpq - stp))
        stpf = stpc;
      else
        stpf = stpq;
      stpf = std::min(stpmax, stpf);
      stpf = std::max(stpmin, stpf);
    }
  } else if (brackt) {
    double theta = three * (fp - fy) / (sty - stp) + dy + dp;
    double s = std::max(std::max(std::fabs(theta), std::fabs(dy)), std::fabs(dp));
    double gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
    if (stp > sty) gamma = -gamma;
    double p = (gamma - dp) + theta;
    double q = ((gamma - dp) + gamma) + dy;
    double r = p / q;
    double stpc = stp + r * (sty - stp);
    stpf = stpc;
  } else if (stp > stx) {
    stpf = stpmax;
  } else {
    stpf = stpmin;
  }
  if (fp > fx) {
    sty = stp;
    fy = fp;
    dy = dp;
  } else {
    if (sgnd < zero) {
      sty = stx;
      fy = fx;
      dy = dx;
    }
    stx = stp;
    fx = fp;
    dx = dp;
  }
  stp = stpf;
}

// ============================================================================
// dpmeps — machine precision.
// ============================================================================
static double dpmeps() {
  return std::numeric_limits<double>::epsilon();
}

// ============================================================================
// dpofa — Cholesky factorization (upper triangle storage, column-major view).
// ============================================================================
static void dpofa(Mat2D& a, int lda, int n, int& info, int r0, int c0) {
  (void)lda;
  for (int j = 1; j <= n; ++j) {
    info = j;
    double s = 0.0;
    int jm1 = j - 1;
    if (jm1 >= 1) {
      for (int k = 1; k <= jm1; ++k) {
        double t = a.at(r0 + k - 1, c0 + j - 1) -
                   ddot1(k - 1, &a.at(r0, c0 + k - 1), 1, &a.at(r0, c0 + j - 1), 1);
        t = t / a.at(r0 + k - 1, c0 + k - 1);
        a.at(r0 + k - 1, c0 + j - 1) = t;
        s = s + t * t;
      }
    }
    s = a.at(r0 + j - 1, c0 + j - 1) - s;
    if (s <= 0.0) return;
    a.at(r0 + j - 1, c0 + j - 1) = std::sqrt(s);
  }
  info = 0;
}

// ============================================================================
// dtrsl — solve triangular system t*x=b (b is a column vector, 1-based).
// ============================================================================
static void dtrsl(Mat2D& t, int ldt, int n, double* b, int job, int& info) {
  (void)ldt;
  for (info = 1; info <= n; ++info) {
    if (t.at(info, info) == 0.0) return;
  }
  info = 0;
  int case_id = 1;
  if (job % 10 != 0) case_id = 2;
  if ((job % 100) / 10 != 0) case_id = case_id + 2;
  if (case_id == 2) {
    b[n - 1] = b[n - 1] / t.at(n, n);
    if (n >= 2) {
      for (int jj = 2; jj <= n; ++jj) {
        int j = n - jj + 1;
        double temp = -b[j + 1 - 1];
        daxpy1(j, temp, &t.at(1, j + 1), 1, b, 1);
        b[j - 1] = b[j - 1] / t.at(j, j);
      }
    }
  } else if (case_id == 3) {
    b[n - 1] = b[n - 1] / t.at(n, n);
    if (n >= 2) {
      for (int jj = 2; jj <= n; ++jj) {
        int j = n - jj + 1;
        b[j - 1] = b[j - 1] - ddot1(jj - 1, &t.at(j + 1, j), 1, &b[j - 1], 1);
        b[j - 1] = b[j - 1] / t.at(j, j);
      }
    }
  } else if (case_id == 4) {
    b[1 - 1] = b[1 - 1] / t.at(1, 1);
    if (n >= 2) {
      for (int j = 2; j <= n; ++j) {
        b[j - 1] = b[j - 1] - ddot1(j - 1, &t.at(1, j), 1, b, 1);
        b[j - 1] = b[j - 1] / t.at(j, j);
      }
    }
  } else {
    b[1 - 1] = b[1 - 1] / t.at(1, 1);
    if (n >= 2) {
      for (int j = 2; j <= n; ++j) {
        double temp = -b[j - 1];
        daxpy1(n - j + 1, temp, &t.at(j, j - 1), 1, &b[j - 1], 1);
        b[j - 1] = b[j - 1] / t.at(j, j);
      }
    }
  }
}
