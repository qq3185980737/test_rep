// ef.cpp — C++ translation of MOPAC 2016 "ef.F90" (2279 lines).
// EF optimizer: eigenvector-following (P-RFO/QA) geometry optimization.
// Subroutines: ef (driver), efsav (restart I/O), efstr (init), formd (step),
// gethes (numerical Hessian), overlp (mode following), prjfc (T/R projection),
// prthes (print), updhes (Powell/BFGS update).
// Conventions: matrices are column-major flat vectors, 1-based indexing
// (index i+(j-1)*nvar = row i, column j), same as common_arrays_C::hesinv.
#include "ef.h"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#ifdef MOPAC_USE_MKL
#include <mkl.h>
#endif
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "ef_C.h"
#include "maps_C.h"
#include "molkst_C.h"
using namespace ef_C;
using namespace common_arrays_C;
using namespace molkst_C;
using chanel_C::iw;
using chanel_C::iw0;
using namespace maps_C;
// External symbols (MOPAC core / BLAS / LAPACK).
extern double ddot(int n, const double* x, int incx, const double* y, int incy);
extern void compfg(const std::vector<double>& xparam, bool int_flag,
                   double& escf, bool fulscf, std::vector<double>& grad,
                   bool lgrad);
extern void symtry();
extern void prttim(double tleft, double& tprt, char& txt);
extern double second(int n);
extern double reada(const std::string& line, int istart);
extern void geout(int iprt);  // must match linker's geout (extern "C" in stub)
extern void mopend(const char* message);
extern void den_in_out(int mode);
extern void prtgra();
extern void write_cell(int iprt);
extern void to_screen(const std::string& text);
extern void rsp(double* a, int n, double* root, double* vect);
extern void dgefa(double* a, int lda, int n, int* ipvt, int& info);
extern void dgedi(double* a, int lda, int n, const int* ipvt, double* det,
                  double* work, int job);
#ifdef MOPAC_USE_MKL
extern "C" void dgemm_(const char* transa, const char* transb, const int* m,
                       const int* n, const int* k, const double* alpha,
                       const double* a, const int* lda, const double* b,
                       const int* ldb, const double* beta, double* c,
                       const int* ldc);
#endif
// Constants from ef.F90 parameter block.
static const double demin = 1.0e-2, gmin = 5.0;
// ---------------------------------------------------------------------------
// helper: C = A*B for n x n column-major flat matrices (beta=0 overwrite).
// Official ef.F90 prjfc calls MKL dgemm; use the same BLAS so the reduction
// order (and hence every bit of f = P*(F*P)) matches the reference exe.
// ---------------------------------------------------------------------------
static void gemm_cpu(int n, const double* a, const double* b, double* c) {
#ifdef MOPAC_USE_MKL
    static int dbg_done = 0;
    if (!dbg_done) {
        dbg_done = 1;
        const char* p = std::getenv("MOPAC_GEMM_DBG");
        if (p) std::fprintf(stderr,
            "DBG gemm n=%d max_threads=%d omp=%s mkl=%s dyn=%s\n",
            n, mkl_get_max_threads(),
            std::getenv("OMP_NUM_THREADS") ? std::getenv("OMP_NUM_THREADS") : "-",
            std::getenv("MKL_NUM_THREADS") ? std::getenv("MKL_NUM_THREADS") : "-",
            std::getenv("MKL_DYNAMIC") ? std::getenv("MKL_DYNAMIC") : "-");
    }
    // The prjfc buffers (hesinv/p/cof) are 1-based-layout vectors: slot 0 is a
    // dead pad, real data is [1..n*n].  Fortran dgemm reads from the address
    // given, so pass data()+1 to skip the pad and align exactly with the
    // reference exe's arrays.  Do NOT pin the thread count: the reference exe
    // runs with MKL_NUM_THREADS=4 and its dgemm (n=114, ~3M FLOPs) goes
    // parallel, so the reduction order must match that same thread count.
    char N = 'N';
    double one = 1.0, zero = 0.0;
    // Reference exe shows NThr:4 in MKL_VERBOSE, but forcing 4 threads here
    // made the trajectory worse (255 steps vs 294), while leaving MKL to its
    // own single-thread path gave 296 steps -- so the reference dgemm is in
    // fact single-threaded for these small matrices.  Keep it unpinned.
    dgemm_(&N, &N, &n, &n, &n, &one, a + 1, &n, b + 1, &n, &zero, c + 1, &n);
#else
        for (int i = 1; i <= n; ++i) c[(j - 1) * n + i] = 0.0;
    for (int j = 1; j <= n; ++j) {
        int cbase = (j - 1) * n;
        for (int k = 1; k <= n; ++k) {
            double akj = b[(j - 1) * n + k];  // b(k,j)
            if (akj == 0.0) continue;
            int abase = (k - 1) * n;
            int i = 1;
            // Manual 4-way unroll: helps MSVC keep 4 independent FMAs in
            // flight for these small (nvar<=~100) matrices.
            for (; i + 3 <= n; i += 4) {
                c[cbase + i] += a[abase + i] * akj;
                c[cbase + i + 1] += a[abase + i + 1] * akj;
                c[cbase + i + 2] += a[abase + i + 2] * akj;
                c[cbase + i + 3] += a[abase + i + 3] * akj;
            }
            for (; i <= n; ++i) c[cbase + i] += a[abase + i] * akj;
        }
    }
#endif  // MOPAC_USE_MKL
}
// ---------------------------------------------------------------------------
// efsav: store/retrieve EF restart data (ef.F90 lines 742-869).
// ---------------------------------------------------------------------------
static void efsav(double tt0, std::vector<double>& hess, double funct,
                  std::vector<double>& grad, std::vector<double>& xparam,
                  std::vector<double>& pmat, int& il, std::vector<double>& bmat,
                  std::vector<int>& ipow, std::vector<double>& oldf,
                  std::vector<double>& d, std::vector<double>& vmode) {
    const int nvar = molkst_C::nvar;
    if (is_PARAM) return;
    const int linear = (nvar * (nvar + 1)) / 2;
    if (ipow[9] == 1 || ipow[9] == 2) {
        // ---- dump ----
        double funct1 = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        if (ipow[9] == 1 && keywrd.find("STEP1") == std::string::npos) {
            std::printf("\n\n          CURRENT VALUE OF GRADIENT NORM =%12.6f\n",
                        funct1);
            if (prt_gradients && keywrd.find(" GRADI") != std::string::npos &&
                mozyme) {
                std::printf("\n\n\n          CURRENT  POINT  AND  DERIVATIVES\n");
                prtgra();
            }
            std::printf("\n          CURRENT VALUE OF GEOMETRY\n");
            geout(iw);
        }
        ipow[8] = nscf;
        std::ofstream of(chanel_C::restart_fn,
                         std::ios::binary | std::ios::trunc);
        of.write(reinterpret_cast<const char*>(&numat), sizeof(int));
        of.write(reinterpret_cast<const char*>(&norbs), sizeof(int));
        for (int i = 1; i <= nvar; ++i)
            of.write(reinterpret_cast<const char*>(&xparam[i]), sizeof(double));
        if (latom != 0) {
            if (keywrd.find(" STEP=") != std::string::npos) {
                of.write(reinterpret_cast<const char*>(&kloop), sizeof(int));
                of.write(reinterpret_cast<const char*>(&rxn_coord),
                         sizeof(double));
                for (int i = 1; i <= kloop; ++i)
                    of.write(reinterpret_cast<const char*>(&profil[i]),
                             sizeof(double));
            } else {
                for (int j = 1; j <= 3; ++j)
                    for (int i = 1; i <= nvar; ++i) {
                        double a = (j <= (int)alparm.size() && i <= (int)alparm[j - 1].size())
                                       ? alparm[j - 1][i - 1] : 0.0;
                        of.write(reinterpret_cast<const char*>(&a),
                                 sizeof(double));
                    }
                of.write(reinterpret_cast<const char*>(&iloop), sizeof(int));
                of.write(reinterpret_cast<const char*>(&x0), sizeof(double));
                of.write(reinterpret_cast<const char*>(&x1), sizeof(double));
                of.write(reinterpret_cast<const char*>(&x2), sizeof(double));
            }
        }
        for (int i = 1; i <= 9; ++i)
            of.write(reinterpret_cast<const char*>(&ipow[i]), sizeof(int));
        of.write(reinterpret_cast<const char*>(&il), sizeof(int));
        of.write(reinterpret_cast<const char*>(&nstep), sizeof(int));
        of.write(reinterpret_cast<const char*>(&funct), sizeof(double));
        of.write(reinterpret_cast<const char*>(&tt0), sizeof(double));
        for (int i = 1; i <= nvar; ++i)
            of.write(reinterpret_cast<const char*>(&grad[i]), sizeof(double));
        for (int j = 1; j <= nvar; ++j)
            for (int i = 1; i <= nvar; ++i)
                of.write(reinterpret_cast<const char*>(&hess[i + (j - 1) * nvar]),
                         sizeof(double));
        for (int j = 1; j <= nvar; ++j)
            for (int i = 1; i <= nvar; ++i)
                of.write(reinterpret_cast<const char*>(&bmat[i + (j - 1) * nvar]),
                         sizeof(double));
        for (int i = 1; i <= nvar; ++i)
            of.write(reinterpret_cast<const char*>(&oldf[i]), sizeof(double));
        for (int i = 1; i <= nvar; ++i)
            of.write(reinterpret_cast<const char*>(&d[i]), sizeof(double));
        for (int i = 1; i <= nvar; ++i)
            of.write(reinterpret_cast<const char*>(&vmode[i]), sizeof(double));
        of.write(reinterpret_cast<const char*>(&ddx), sizeof(double));
        of.write(reinterpret_cast<const char*>(&ef_mode), sizeof(int));
        of.write(reinterpret_cast<const char*>(&nstep), sizeof(int));
        of.write(reinterpret_cast<const char*>(&negreq), sizeof(int));
        for (int i = 1; i <= linear; ++i)
            of.write(reinterpret_cast<const char*>(&pmat[i]), sizeof(double));
        den_in_out(1);
        if (keywrd.find("STEP1") != std::string::npos) return;
        // of closed on destruction
        return;
    }
    // ---- restore ----
    std::ifstream inf(chanel_C::restart_fn, std::ios::binary);
    if (!inf.good()) {
        mopend("NO RESTART FILE EXISTS!");
        return;
    }
    int old_numat = 0, old_norbs = 0;
    bool io_ok = inf.read(reinterpret_cast<char*>(&old_numat), sizeof(int))
                     .good() &&
                 inf.read(reinterpret_cast<char*>(&old_norbs), sizeof(int))
                     .good();
    for (int i = 1; i <= nvar && io_ok; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&xparam[i]), sizeof(double))
                    .good();
    if (!io_ok) {
        mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
        return;
    }
    if (norbs != old_norbs || numat != old_numat) {
        mopend("Restart file read in does not match current data set");
        return;
    }
    if (latom != 0) {
        if (keywrd.find(" STEP=") != std::string::npos) {
            io_ok = inf.read(reinterpret_cast<char*>(&kloop), sizeof(int))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&rxn_coord), sizeof(double))
                        .good();
            for (int i = 1; i <= kloop && io_ok; ++i)
                io_ok = inf.read(reinterpret_cast<char*>(&profil[i]),
                                 sizeof(double)).good();
        } else {
            for (int j = 1; j <= 3 && io_ok; ++j)
                for (int i = 1; i <= nvar && io_ok; ++i) {
                    double a = 0.0;
                    io_ok = inf.read(reinterpret_cast<char*>(&a), sizeof(double))
                                .good();
                    if (j <= (int)alparm.size() &&
                        i <= (int)alparm[j - 1].size())
                        alparm[j - 1][i - 1] = a;
                }
            io_ok = inf.read(reinterpret_cast<char*>(&iloop), sizeof(int))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&x0), sizeof(double))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&x1), sizeof(double))
                        .good() &&
                    inf.read(reinterpret_cast<char*>(&x2), sizeof(double))
                        .good();
        }
        if (!io_ok) {
            mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
            return;
        }
    }
    for (int i = 1; i <= 9; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&ipow[i]), sizeof(int)).good() &&
                io_ok;
    io_ok = inf.read(reinterpret_cast<char*>(&il), sizeof(int)).good() &&
            inf.read(reinterpret_cast<char*>(&nstep), sizeof(int)).good() &&
            inf.read(reinterpret_cast<char*>(&funct), sizeof(double)).good() &&
            inf.read(reinterpret_cast<char*>(&tt0), sizeof(double)).good();
    if (!io_ok) {
        mopend("RESTART FILE EXISTS, BUT IS CORRUPT");
        return;
    }
    nscf = ipow[8];
    int k = (int)(tt0 / 1000000.0);
    tt0 = tt0 - k * 1000000.0;
    std::printf("\n          TOTAL TIME USED SO FAR:%13.2f SECONDS\n", tt0);
    if (std::abs(funct) > 1.0e-20)
        std::printf("                         FUNCTION:%17.6f\n", funct);
    for (int i = 1; i <= nvar; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&grad[i]), sizeof(double))
                    .good() && io_ok;
    for (int j = 1; j <= nvar; ++j)
        for (int i = 1; i <= nvar; ++i)
            io_ok = inf
                        .read(reinterpret_cast<char*>(
                                  &hess[i + (j - 1) * nvar]),
                              sizeof(double))
                        .good() && io_ok;
    for (int j = 1; j <= nvar; ++j)
        for (int i = 1; i <= nvar; ++i)
            io_ok = inf
                        .read(reinterpret_cast<char*>(
                                  &bmat[i + (j - 1) * nvar]),
                              sizeof(double))
                        .good() && io_ok;
    for (int i = 1; i <= nvar; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&oldf[i]), sizeof(double))
                    .good() && io_ok;
    for (int i = 1; i <= nvar; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&d[i]), sizeof(double))
                    .good() && io_ok;
    for (int i = 1; i <= nvar; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&vmode[i]), sizeof(double))
                    .good() && io_ok;
    io_ok = inf.read(reinterpret_cast<char*>(&ddx), sizeof(double)).good() &&
            inf.read(reinterpret_cast<char*>(&ef_mode), sizeof(int)).good() &&
            inf.read(reinterpret_cast<char*>(&nstep), sizeof(int)).good() &&
            inf.read(reinterpret_cast<char*>(&negreq), sizeof(int)).good();
    for (int i = 1; i <= linear && io_ok; ++i)
        io_ok = inf.read(reinterpret_cast<char*>(&pmat[i]), sizeof(double))
                    .good() && io_ok;
    if (keywrd.find("STEP1") != std::string::npos) return;
    if (io_ok == false) mopend("Restart file is currupt");
}
// ---------------------------------------------------------------------------
// efstr: initialization (ef.F90 lines 870-1094).
// ---------------------------------------------------------------------------
static void efstr(std::vector<double>& xparam, double& funct, int& ihess,
                  int& ntime, int& iloop, int& igthes, int& mxstep,
                  int& ireclc, int& iupd, double& dmax, double& ddmax,
                  double& ddmin, double& tol2, double& time1, double& time2,
                  int& nvar, bool& scf1, bool& lupd, int& ldump, bool& rrscal,
                  bool& donr, std::vector<double>& hess,
                  std::vector<double>& bmat, std::vector<double>& pmat,
                  std::vector<double>& grad, std::vector<double>& oldf,
                  std::vector<double>& d, std::vector<double>& vmode) {
    static int icalcn = 0;
    bool restrt;
    int i, ip, its, j, k, mtmp;
    double tt0;
    std::vector<int> ipow(10, 0);
    nvar = std::abs(nvar);
    ldump = 0;
    lupd = (keywrd.find(" NOUPD") == std::string::npos);
    restrt = (keywrd.find(" RESTART") != std::string::npos);
    scf1 = (keywrd.find(" 1SCF") != std::string::npos);
    nstep = 0;
    ihess = 0;
    ntime = 0;
    iloop = 1;
    ef_mode = 0;
    igthes = 0;
    iupd = 2;
    negreq = 0;
    rmin = 0.0;
    rmax = 1.0e3;
    dmax = 0.2;
    ddmax = 0.5;
    if (icalcn == numcal) {
        ddmax = 0.1;
        ireclc = 10;
    } else {
        icalcn = numcal;
        ireclc = 999999;
    }
    limscf = true;
    its = (int)keywrd.find(" TS ");
    if (its != (int)std::string::npos) {
        limscf = false;
        ef_mode = 1;
        igthes = 1;
        iupd = 1;
        negreq = 1;
        rmin = 0.0;
        rmax = 4.0;
        omin = 0.8;
        dmax = 0.1;
        ddmax = 0.3;
    }
    rrscal = false;
    i = (int)keywrd.find(" RSCAL");
    if (i != (int)std::string::npos) rrscal = true;
    donr = true;
    i = (int)keywrd.find(" NONR");
    if (i != (int)std::string::npos) donr = false;
    iprnt = 3;
    ip = (int)keywrd.find(" PRNT=");
    if (ip != (int)std::string::npos) {
        iprnt = (int)std::lround(reada(keywrd, ip));
        if (iprnt > 5) iprnt = 5;
        if (iprnt < 0) iprnt = 0;
    }
    mxstep = 2000;
    i = (int)keywrd.find(" CYCLES=");
    if (i != (int)std::string::npos) {
        mxstep = (int)std::lround(reada(keywrd, i));
        if (mxstep == 0 && ip == (int)std::string::npos) iprnt = 3;
    }
    i = (int)keywrd.find(" RECALC=");
    if (i != (int)std::string::npos)
        ireclc = (int)std::lround(reada(keywrd, i));
    i = (int)keywrd.find(" IUPD=");
    if (i != (int)std::string::npos) iupd = (int)std::lround(reada(keywrd, i));
    i = (int)keywrd.find(" MODE=");
    if (i != (int)std::string::npos)
        ef_mode = (int)std::lround(reada(keywrd, i));
    ddmin = 1.0e-4;
    i = (int)keywrd.find(" DDMIN=");
    if (i != (int)std::string::npos) ddmin = reada(keywrd, i);
    i = (int)keywrd.find(" DMAX=");
    if (i != (int)std::string::npos) dmax = reada(keywrd, i);
    i = (int)keywrd.find(" DDMAX=");
    if (i != (int)std::string::npos) ddmax = reada(keywrd, i);
    tol2 = 1.0;
    if (id != 0) tol2 = id * 2.0 - 1.0;
    if (keywrd.find(" PREC") != std::string::npos) tol2 *= 5.0e-2;
    i = (int)keywrd.find(" GNORM=");
    if (i != (int)std::string::npos) tol2 = reada(keywrd, i);
    if (keywrd.find(" LET") == std::string::npos && tol2 < 0.01) {
        std::printf("\n  GNORM HAS BEEN SET TOO LOW, RESET TO 0.01\n"
                    "  SPECIFY LET AS KEYWORD TO ALLOW GNORM LESS THAN 0.01\n");
        tol2 = 0.01;
    }
    i = (int)keywrd.find(" HESS=");
    if (i != (int)std::string::npos) igthes = (int)std::lround(reada(keywrd, i));
    i = (int)keywrd.find(" RMIN=");
    if (i != (int)std::string::npos) rmin = reada(keywrd, i);
    i = (int)keywrd.find(" RMAX=");
    if (i != (int)std::string::npos) rmax = reada(keywrd, i);
    i = (int)keywrd.find(" OMIN=");
    if (i != (int)std::string::npos) omin = reada(keywrd, i);
    time1 = time0;
    time2 = time1;
    if ((its != (int)std::string::npos) && (iupd == 2)) {
        std::printf(" TS SEARCH AND BFGS UPDATE WILL NOT WORK\n");
        mopend("TS SEARCH AND BFGS UPDATE WILL NOT WORK");
        return;
    }
    if ((its != (int)std::string::npos) && (igthes == 0)) {
        std::printf(" TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN\n");
        mopend("TS SEARCH REQUIRE BETTER THAN DIAGONAL HESSIAN");
        return;
    }
    if ((igthes < 0) || (igthes > 3)) {
        std::printf(" UNRECOGNIZED HESS OPTION %d\n", igthes);
        mopend("UNRECOGNIZED HESS OPTION");
        return;
    }
    if ((omin < 0.0) || (omin > 1.0)) {
        std::printf(" OMIN MUST BE BETWEEN 0 AND 1 %g\n", omin);
        mopend("OMIN MUST BE BETWEEN 0 AND 1");
        return;
    }
    nstep = 0;
    if (restrt) {
        ipow[9] = 0;
        mtmp = ef_mode;
        tt0 = 0.0;
        j = 0;
        efsav(tt0, hess, funct, grad, xparam, pmat, i, bmat, ipow, oldf, d,
              vmode);
        if (moperr) return;
        ef_mode = mtmp;
        j = ipow[2];
        k = (int)(tt0 / 1000000.0);
        time0 = time0 - tt0 + k * 1000000.0;
        iloop = i;
        if (i > 0) {
            igthes = 4;
            nstep = j;
            std::printf("          RESTARTING HESSIAN AT POINT%4d\n", iloop);
            if (nstep != 0)
                std::printf("          IN OPTIMIZATION STEP%4d\n", nstep);
        } else {
            nstep = j;
            std::printf("\n          RESTARTING OPTIMIZATION AT STEP%6d\n",
                        nstep);
        }
    }
    mxstep = mxstep + nstep;
}
// ---------------------------------------------------------------------------
// overlp: mode following overlap (ef.F90 lines 1791-1875).
// ---------------------------------------------------------------------------
static void overlp(double dmax, double ddmin, int& newmod, int nvar,
                   bool& lorjk, const std::vector<double>& u,
                   std::vector<double>& vmode) {
    static int icalcn = 0, it = 0;
    double ovlp, tovlp;
    if (icalcn != numcal) {
        icalcn = numcal;
        if (ef_mode > nvar) {
            std::printf("ERROR!! MODE IS LARGER THAN NVAR\n");
            mopend("ERROR!! MODE IS LARGER THAN NVAR");
            return;
        } else {
            it = ef_mode;
            if (iprnt >= 1)
                std::printf(" HESSIAN MODE FOLLOWING SWITCHED ON\n"
                            "     FOLLOWING MODE %3d\n",
                            ef_mode);
        }
    } else {
        it = 1;
        lorjk = false;
        tovlp = 0.0;
        for (int k = 1; k <= nvar; ++k)
            tovlp += u[k + (it - 1) * nvar] * vmode[k];
        tovlp = std::abs(tovlp);
        for (int i = 2; i <= nvar; ++i) {
            ovlp = 0.0;
            for (int k = 1; k <= nvar; ++k)
                ovlp += u[k + (i - 1) * nvar] * vmode[k];
            ovlp = std::abs(ovlp);
            if (ovlp > tovlp) {
                tovlp = ovlp;
                it = i;
            }
        }
        if (iprnt >= 1)
            std::printf("OVERLAP OF CURRENT MODE%3d WITH PREVIOUS MODE IS %6.3f\n",
                        it, tovlp);
        if (tovlp < omin) {
            if (dmax > ddmin) {
                lorjk = true;
                if (iprnt >= 1)
                    std::printf("OVERLAP LESS THAN OMIN%6.3f REJECTING PREVIOUS STEP\n",
                                omin);
                return;
            } else if (iprnt >= 1) {
                std::printf("OVERLAP LESS THAN OMIN%6.3f BUT TRUST RADIUS%6.3f "
                            "IS LESS THAN DDMIN%6.3f\n ACCEPTING STEP\n",
                            omin, dmax, ddmin);
            }
        }
    }
    for (int i = 1; i <= nvar; ++i) vmode[i] = u[i + (it - 1) * nvar];
    newmod = it;
}
// ---------------------------------------------------------------------------
// prthes: print Hessian (ef.F90 lines 2119-2174).
// ---------------------------------------------------------------------------
static void prthes(const std::vector<double>& eigval, int nvar,
                   const std::vector<double>& hess,
                   const std::vector<double>& u) {
    if (iprnt >= 4) {
        std::printf("\n \n              HESSIAN MATRIX\n");
        int low = 1, nup = 8;
        for (;;) {
            nup = std::min(nup, nvar);
            std::printf("\n");
            for (int i = low; i <= nup; ++i) std::printf("%9d", i);
            std::printf("\n");
            for (int i = 1; i <= nvar; ++i) {
                std::printf("%3d", i);
                for (int j = low; j <= nup; ++j)
                    std::printf("%9.1f", hess[i + (j - 1) * nvar]);
                std::printf("\n");
            }
            nup += 8;
            low += 8;
            if (low > nvar) break;
        }
    }
    std::printf("\n \n              HESSIAN EIGENVALUES AND -VECTORS\n");
    int low = 1, nup = 8;
    for (;;) {
        nup = std::min(nup, nvar);
        std::printf("\n");
        for (int i = low; i <= nup; ++i) std::printf("%9d", i);
        std::printf("\n");
        for (int i = low; i <= nup; ++i) std::printf("%9.1f", eigval[i]);
        std::printf("\n");
        for (int i = 1; i <= nvar; ++i) {
            std::printf("%3d", i);
            for (int j = low; j <= nup; ++j)
                std::printf("%9.4f", u[i + (j - 1) * nvar]);
            std::printf("\n");
        }
        nup += 8;
        low += 8;
        if (low > nvar) break;
    }
}
// ---------------------------------------------------------------------------
// updhes: Powell/BFGS Hessian update (ef.F90 lines 2175-2303).
// ---------------------------------------------------------------------------
void updhes(std::vector<double>& svec, std::vector<double>& tvec,
            const std::vector<double>& grad, int nvar, int iupd,
            std::vector<double>& hess, const std::vector<double>& oldf,
            const std::vector<double>& d) {
    static int icalcn = 0;
    if (icalcn != numcal) {
        icalcn = numcal;
        if (iprnt >= 2) {
            if (iupd == 0) std::printf("\n HESSIAN IS NOT BEING UPDATED\n");
            if (iupd == 1)
                std::printf("\n HESSIAN IS BEING UPDATED USING THE POWELL UPDATE\n");
            if (iupd == 2)
                std::printf("\n HESSIAN IS BEING UPDATED USING THE BFGS UPDATE\n");
        }
    }
    if (iupd == 0) return;
    for (int i = 1; i <= nvar; ++i) tvec[i] = 0.0;
    for (int j = 1; j <= nvar; ++j)
        for (int i = 1; i <= nvar; ++i)
            tvec[i] += hess[i + (j - 1) * nvar] * d[j];
    double dds = 0.0, ddtd = 0.0, temp;
    if (iupd == 1) {
        for (int i = 1; i <= nvar; ++i) {
            tvec[i] = grad[i] - oldf[i] - tvec[i];
            svec[i] = grad[i] - oldf[i];
        }
        dds = ddx * ddx;
        ddtd = 0.0;
        for (int i = 1; i <= nvar; ++i) ddtd += tvec[i] * d[i];
        if (std::abs(dds) < 1.0e-20) {
            std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
            mopend("in UPDHES");
            return;
        }
        ddtd = ddtd / dds;
        for (int i = 2; i <= nvar; ++i)
            for (int j = 1; j < i; ++j) {
                temp = tvec[i] * d[j] + d[i] * tvec[j] - d[i] * ddtd * d[j];
                hess[i + (j - 1) * nvar] += temp / dds;
                hess[j + (i - 1) * nvar] = hess[i + (j - 1) * nvar];
            }
        for (int i = 1; i <= nvar; ++i) {
            temp = d[i] * (2.0 * tvec[i] - d[i] * ddtd);
            hess[i + (i - 1) * nvar] += temp / dds;
        }
        return;
    }
    if (iupd != 2) return;
    for (int i = 1; i <= nvar; ++i) svec[i] = grad[i] - oldf[i];
    dds = 0.0;
    ddtd = 0.0;
    for (int i = 1; i <= nvar; ++i) {
        dds += svec[i] * d[i];
        ddtd += d[i] * tvec[i];
    }
    if (std::abs(dds) < 1.0e-20 || std::abs(ddtd) < 1.0e-20) {
        std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
        mopend("in UPDHES");
        return;
    }
    for (int i = 2; i <= nvar; ++i)
        for (int j = 1; j < i; ++j) {
            temp = (svec[i] * svec[j]) / dds - (tvec[i] * tvec[j]) / ddtd;
            hess[i + (j - 1) * nvar] += temp;
            hess[j + (i - 1) * nvar] = hess[i + (j - 1) * nvar];
        }
    for (int i = 1; i <= nvar; ++i) {
        temp = (svec[i] * svec[i]) / dds - (tvec[i] * tvec[i]) / ddtd;
        hess[i + (i - 1) * nvar] += temp;
    }
}
// ---------------------------------------------------------------------------
// gethes: numerical Hessian (ef.F90 lines 1535-1790).
// ---------------------------------------------------------------------------
static void gethes(std::vector<double>& xparam, int igthes, int& iloop,
                   std::vector<double>& hess, std::vector<double>& pmat,
                   std::vector<double>& bmat, std::vector<double>& grad,
                   std::vector<std::vector<double>>& geo,
                   const std::vector<std::vector<int>>& loc,
                   std::vector<double>& oldf, std::vector<double>& d,
                   std::vector<double>& vmode, double& funct0) {
    const int nvar = molkst_C::nvar;
    const double dghsa = 500.0, dghsd = 200.0, dghss = 1000.0;
    const double two = 2.0, xinc = 1.0e-3, zzero = 0.0;
    std::vector<double> gnext1(nvar + 1, 0.0), gmin1(nvar + 1, 0.0);
    std::vector<int> ipow(10, 0);
    int i, j, k, l, mtmp, nxxx, ij, percent;
    double dg1, dg2, dg3, dummy, fdmy, funct = 0.0, funct1, sum, test = 0.0, tdm,
        time1, time2, tstep, tstore, tt0;
    bool lpacifier;
    (void)pmat;
    (void)bmat;
    (void)d;
    (void)vmode;
    if (igthes == 0) {
        std::printf("\n          DIAGONAL MATRIX USED AS START HESSIAN\n");
        for (i = 1; i <= nvar; ++i)
            for (j = 1; j <= nvar; ++j) hess[i + (j - 1) * nvar] = zzero;
        if (id != 0) {
            for (i = 1; i <= nvar; ++i)
                oldf[i] = xparam[i] - (grad[i] < 0.0 ? -0.001 : 0.001);
            compfg(oldf, true, funct1, true, gnext1, true);
            for (i = 1; i <= nvar; ++i)
                hess[i + (i - 1) * nvar] = std::max(
                    100.0, (grad[i] - gnext1[i]) / (xparam[i] - oldf[i]));
            if (funct1 < funct0) {
                funct0 = funct1;
                for (i = 1; i <= nvar; ++i) {
                    xparam[i] = oldf[i];
                    grad[i] = gnext1[i];
                }
            }
        } else {
            sum = 200.0;
            if (keywrd.find(" XYZ") != std::string::npos) {
                dg1 = sum;
                dg2 = sum;
                dg3 = sum;
            } else {
                dg1 = dghss;
                dg2 = dghsa;
                dg3 = dghsd;
            }
            ij = 1;
            for (j = 1; j <= natoms; ++j) {
                for (i = 1; i <= 3; ++i) {
                    if (loc[2][ij] == i && loc[1][ij] == j) {
                        if (i == 1) hess[ij + (ij - 1) * nvar] = dg1;
                        if (i == 2) hess[ij + (ij - 1) * nvar] = dg2;
                        if (i == 3) hess[ij + (ij - 1) * nvar] = dg3;
                        ++ij;
                        if (ij > nvar) break;
                    }
                }
                if (ij > nvar) break;
            }
            ij = ij - 1;
            if (ij != nvar)
                std::printf("ERROR IN IGTHES=0,IJ,NVAR %d %d\n", ij, nvar);
        }
    }
    if (igthes == 2) {
        std::printf("\n          HESSIAN READ FROM DISK\n");
        ipow[9] = 0;
        nxxx = nalpha;
        nalpha = 0;
        mtmp = ef_mode;
        tdm = 0.0;
        fdmy = 0.0;
        int iidum = 0;
        efsav(tdm, hess, fdmy, gnext1, gmin1, pmat, iidum, bmat, ipow, oldf, d,
              vmode);
        if (moperr) return;
        nalpha = nxxx;
        ef_mode = mtmp;
        nstep = 0;
    }
    if ((igthes == 1) || (igthes == 3) || (igthes == 4)) {
        if (igthes == 1) std::printf("\n          HESSIAN CALCULATED NUMERICALLY\n");
        if (igthes == 3)
            std::printf("\n          HESSIAN CALCULATED DOUBLE NUMERICALLY\n");
        if (iprnt >= 5) {
            std::printf("%3d", 0);
            for (int jj = 1; jj <= nvar; ++jj) std::printf("%9.4f", grad[jj]);
            std::printf("\n");
        }
        time1 = second(1);
        tstore = time1;
        percent = (100 * (iloop - 1)) / nvar;
        lpacifier = false;
        for (i = iloop; i <= nvar; ++i) {
            xparam[i] = xparam[i] + xinc;
            compfg(xparam, true, dummy, true, gnext1, true);
            if (iprnt >= 5) {
                std::printf("%3d", i);
                for (int jj = 1; jj <= nvar; ++jj)
                    std::printf("%9.4f", gnext1[jj]);
                std::printf("\n");
            }
            xparam[i] = xparam[i] - xinc;
            if (igthes == 3) {
                xparam[i] = xparam[i] - xinc;
                gnorm = 0.0;
                compfg(xparam, true, dummy, true, gmin1, true);
                if (iprnt >= 5) {
                    std::printf("%3d", -i);
                    for (int jj = 1; jj <= nvar; ++jj)
                        std::printf("%9.4f", gmin1[jj]);
                    std::printf("\n");
                }
                xparam[i] = xparam[i] + xinc;
                for (j = 1; j <= nvar; ++j)
                    hess[i + (j - 1) * nvar] = (gnext1[j] - gmin1[j]) / (xinc + xinc);
            } else {
                for (j = 1; j <= nvar; ++j)
                    hess[i + (j - 1) * nvar] = (gnext1[j] - grad[j]) / xinc;
            }
            time2 = second(1);
            tstep = time2 - time1;
            tleft = tleft - tstep;
            time1 = time2;
            j = (100 * i) / nvar;
            test = test + tstep;
            if (j / 5 > percent / 5 || (j == 1 && percent != 1) || i == 1) {
                percent = j;
                j = i - iloop + 1;
                if (nvar * test / (j * 3600) > 0.1 || lpacifier) {
                    lpacifier = true;
                    sum = (nvar - i) * (test / (j * 3600));
                    char buf[128];
                    std::snprintf(buf, sizeof(buf),
                                  "    Hessian%3d%% complete.  Estimated "
                                  "remaining time required:%7.2f hours",
                                  percent, sum);
                    line = buf;
                    to_screen(line);
                }
            }
            if (nvar * test / (j * 3600) > 0.1 || lpacifier) {
                char buf[128];
                std::snprintf(buf, sizeof(buf), "%5d of%4d steps completed", i,
                              nvar);
                line = buf;
                to_screen(line);
            }
            if (tleft < tstep * two) {
                std::printf(" NOT ENOUGH TIME TO COMPLETE HESSIAN\n");
                std::printf(" STOPPING IN HESSIAN AT COORDINATE: %d\n", i);
                ipow[9] = 1;
                ipow[2] = nstep;
                tt0 = second(1) - time0;
                j = i;
                efsav(tt0, hess, funct, grad, xparam, pmat, j, bmat, ipow, oldf,
                      d, vmode);
                tleft = -0.1;
                mopend("NOT ENOUGH TIME TO COMPLETE HESSIAN");
                return;
            }
        }
        k = loc[1][nvar];
        l = loc[2][nvar];
        geo[l][k] = xparam[nvar];
        if (ndep != 0) symtry();
        time2 = second(1);
        tstep = time2 - tstore;
        tleft = tleft + tstep;
    }
    for (i = 1; i <= nvar; ++i)
        for (j = 1; j < i; ++j) {
            hess[i + (j - 1) * nvar] =
                (hess[i + (j - 1) * nvar] + hess[j + (i - 1) * nvar]) / two;
            hess[j + (i - 1) * nvar] = hess[i + (j - 1) * nvar];
        }
}
// ---------------------------------------------------------------------------
// prjfc: remove translations/rotations for Cartesian optimization
// (ef.F90 lines 1876-2118).
// ---------------------------------------------------------------------------
static void prjfc(std::vector<double>& f, const std::vector<double>& xparam,
                  int nvar, std::vector<double>& cof, std::vector<double>& p,
                  const std::vector<double>& atmass, std::vector<double>& x,
                  std::vector<double>& rm, std::vector<double>& dx,
                  std::vector<std::vector<double>>& coord) {
    const double zero = 0.0, one = 1.0, cut8 = 1.0e-8;
    // Levi-Civita-like tensor, Fortran data order (ic fastest: ib, then ia).
    static const double tens[27] = {
        0, 0, 0, 0, 0, -1, 0, 1, 0,      // ic=1
        0, 0, 1, 0, 0, 0, -1, 0, 0,      // ic=2
        0, -1, 0, 1, 0, 0, 0, 0, 0       // ic=3
    };
    int natm_loc = nvar / 3;
    int nc1 = nvar;
    int ij = 1;
    for (int i = 1; i <= natm_loc; ++i) {
        coord[0][i] = xparam[ij];
        coord[1][i] = xparam[ij + 1];
        coord[2][i] = xparam[ij + 2];
        ij = ij + 3;
    }
    int l = 0;
    for (int i = 1; i <= natm_loc; ++i) {
        double tmp = one / std::sqrt(atmass[i]);
        for (int j = 1; j <= 3; ++j) rm[++l] = tmp;
    }
    for (int i = 1; i <= nc1; ++i) dx[i] = zero;
    double totm = zero;
    double cmass[4] = {0.0, 0.0, 0.0, 0.0};
    for (int i = 1; i <= natm_loc; ++i) {
        totm = totm + atmass[i];
        for (int j = 1; j <= 3; ++j) cmass[j] += atmass[i] * coord[j - 1][i];
    }
    for (int j = 1; j <= 3; ++j) cmass[j] = cmass[j] / totm;
    l = 0;
    for (int i = 1; i <= natm_loc; ++i)
        for (int j = 1; j <= 3; ++j) {
            double tmp = std::sqrt(atmass[i]);
            x[++l] = tmp * (coord[j - 1][i] - cmass[j]);
        }
    // Inertia tensor (row-major view rot[i][j] = rot(i,j)).
    // IMPORTANT: Fortran `rot = rot + b + c` is left-assoc ((rot+b)+c);
    // C++ `rot += b + c` computes (b+c) first. Split diagonals into two
    // += so the rounding path matches the reference bit-for-bit.
    double rot[4][4] = {};
    for (int i = 1; i <= natm_loc; ++i) {
        l = 3 * (i - 1) + 1;
        rot[1][1] += x[l + 1] * x[l + 1];
        rot[1][1] += x[l + 2] * x[l + 2];
        rot[1][2] -= x[l] * x[l + 1];
        rot[1][3] -= x[l] * x[l + 2];
        rot[2][2] += x[l] * x[l];
        rot[2][2] += x[l + 2] * x[l + 2];
        rot[2][3] -= x[l + 1] * x[l + 2];
        rot[3][3] += x[l] * x[l];
        rot[3][3] += x[l + 1] * x[l + 1];
    }
    rot[2][1] = rot[1][2];
    rot[3][1] = rot[1][3];
    rot[3][2] = rot[2][3];
    if (std::getenv("PRJFC_PROBE")) {
        static bool rot_dumped = false;
        if (!rot_dumped) {
            rot_dumped = true;
            FILE* fp = std::fopen("prjfc_rot_our.txt", "w");
            if (fp) {
                std::fprintf(fp, "%.17e %.17e %.17e\n", cmass[1], cmass[2],
                             cmass[3]);
                std::fprintf(fp, "%.17e\n", totm);
                for (int ii = 1; ii <= 3; ++ii)
                    std::fprintf(fp, " % .17e % .17e % .17e\n", rot[1][ii],
                                 rot[2][ii], rot[3][ii]);
                for (int ii = 1; ii <= nc1; ++ii)
                    std::fprintf(fp, "%.17e\n", x[ii]);
                std::fclose(fp);
            }
        }
    }
    double chk = rot[1][1] * rot[2][2] * rot[3][3];
    int iscr[6] = {};
    double det2[2] = {};
    double scr[10] = {};
    if (std::abs(chk) > cut8) {
        int info = 0;
        double rot1[9];
        for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= 3; ++j)
                rot1[(j - 1) * 3 + (i - 1)] = rot[i][j];
        dgefa(rot1, 3, 3, iscr, info);
        if (info != 0) mopend("Error in EF");
        if (info != 0) return;
        dgedi(rot1, 3, 3, iscr, det2, scr, 1);
        for (int i = 1; i <= 3; ++i)
            for (int j = 1; j <= 3; ++j)
                rot[i][j] = rot1[(j - 1) * 3 + (i - 1)];
    } else if (std::abs(rot[1][1]) > cut8) {
        if (std::abs(rot[2][2]) > cut8) {
            double det = rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1];
            double trp = rot[1][1];
            rot[1][1] = rot[2][2] / det;
            rot[2][2] = trp / det;
            rot[1][2] = -rot[1][2] / det;
            rot[2][1] = -rot[2][1] / det;
        } else if (std::abs(rot[3][3]) > cut8) {
            double det = rot[1][1] * rot[3][3] - rot[1][3] * rot[3][1];
            double trp = rot[1][1];
            rot[1][1] = rot[3][3] / det;
            rot[3][3] = trp / det;
            rot[1][3] = -rot[1][3] / det;
            rot[3][1] = -rot[3][1] / det;
        } else {
            rot[1][1] = one / rot[1][1];
        }
    } else if (std::abs(rot[2][2]) > cut8) {
        if (std::abs(rot[3][3]) > cut8) {
            double det = rot[3][3] * rot[2][2] - rot[3][2] * rot[2][3];
            double trp = rot[3][3];
            rot[3][3] = rot[2][2] / det;
            rot[2][2] = trp / det;
            rot[3][2] = -rot[3][2] / det;
            rot[2][3] = -rot[2][3] / det;
        } else {
            rot[2][2] = one / rot[2][2];
        }
    } else if (std::abs(rot[3][3]) > cut8) {
        rot[3][3] = one / rot[3][3];
    } else {
        std::printf("EVERY DIAGONAL ELEMENTS ARE ZERO ? %20.10f %20.10f %20.10f\n",
                    rot[1][1], rot[2][2], rot[3][3]);
        return;
    }
    // P matrix.
    // Nonzero (ia,ib) pattern of tens(ia,ib,ic) with signs, in the EXACT
    // enumeration order of the reference loop
    //   do ia = 1,3 ; do ib = 1,3 ; if (tens(ia,ib,ic)/=0) ...
    // (ia outer, ib inner, first-hit order).
    //   ic=1: tens(2,3,1)=+1 then tens(3,2,1)=-1
    //   ic=2: tens(1,3,2)=-1 then tens(3,1,2)=+1
    //   ic=3: tens(1,2,3)=+1 then tens(2,1,3)=-1
    static const int nz_ia[3][2] = {{2, 3}, {1, 3}, {1, 2}};
    static const int nz_ib[3][2] = {{3, 2}, {3, 1}, {2, 1}};
    static const double nz_s[3][2] = {{1.0, -1.0}, {-1.0, 1.0}, {1.0, -1.0}};
    for (int ip = 1; ip <= natm_loc; ++ip) {
        int indx = 3 * (ip - 1);
        for (int jp = 1; jp <= ip; ++jp) {
            int jndx = 3 * (jp - 1);
            for (int ic = 1; ic <= 3; ++ic) {
                int jend = 3;
                if (jp == ip) jend = ic;
                for (int jc = 1; jc <= jend; ++jc) {
                    double sum = 0.0;
                    for (int e1 = 0; e1 < 2; ++e1)
                        for (int e2 = 0; e2 < 2; ++e2)
                            sum += nz_s[ic - 1][e1] * nz_s[jc - 1][e2] *
                                   rot[nz_ia[ic - 1][e1]][nz_ia[jc - 1][e2]] *
                                   x[indx + nz_ib[ic - 1][e1]] *
                                   x[jndx + nz_ib[jc - 1][e2]];
                    int ii = indx + ic;
                    int jj = jndx + jc;
                    p[ii + (jj - 1) * nc1] = sum + dx[ii] * dx[jj];
                    if (ic == jc)
                        p[ii + (jj - 1) * nc1] +=
                            one / (rm[ii] * rm[jj] * totm);
                }
            }
        }
    }
    // delta(i,j) - p(i,j).
    for (int i = 1; i <= nc1; ++i)
        for (int j = 1; j <= i; ++j) {
            p[i + (j - 1) * nc1] = -p[i + (j - 1) * nc1];
            if (i == j) p[i + (j - 1) * nc1] = one + p[i + (j - 1) * nc1];
        }
    // Zero out tiny entries.
    for (int i = 1; i <= nc1; ++i)
        for (int j = 1; j <= i; ++j) {
            if (std::abs(p[i + (j - 1) * nc1]) < cut8)
                p[i + (j - 1) * nc1] = zero;
            p[j + (i - 1) * nc1] = p[i + (j - 1) * nc1];
        }
    // cof = f*p;  f = p*cof   (CPU dgemm equivalents).
    if (std::getenv("PRJFC_PROBE")) {
        // First invocation only: dump p (dgemm input) and f (hessc output).
        static bool dumped = false;
        if (!dumped) {
            dumped = true;
            std::printf("PRJFC cyc nvar=%d\n", nc1);
            std::printf("PRJFC p[1..8] =");
            for (int k = 1; k <= 8; ++k)
                std::printf(" %20.14f", p[k + (0) * nc1]);
            std::printf("\nPRJFC p[%d]= (i=%d,j=1)\n", nc1, nc1);
            std::printf("PRJFC p[i=%d]= %20.14f %20.14f %20.14f\n", nc1,
                        p[nc1 + 0 * nc1], p[nc1 + 1 * nc1], p[nc1 + 2 * nc1]);
            FILE* pf = std::fopen("prjfc_p_our.txt", "w");
            if (pf) {
                for (int ii = 1; ii <= nc1; ++ii)
                    for (int jj = 1; jj <= ii; ++jj)
                        std::fprintf(pf, "%.17e\n", p[jj + (ii - 1) * nc1]);
                std::fclose(pf);
            }
            FILE* hf = std::fopen("prjfc_hes_our.txt", "w");
            if (hf) {
                for (int ii = 1; ii <= nc1; ++ii)
                    for (int jj = 1; jj <= ii; ++jj)
                        std::fprintf(hf, "%.17e\n", f[jj + (ii - 1) * nc1]);
                std::fclose(hf);
            }
            std::printf("PRJFC f[1..8] =");
            for (int k = 1; k <= 8; ++k)
                std::printf(" %20.14f", f[k]);
            std::printf("\n");
        }
    }
    gemm_cpu(nc1, f.data(), p.data(), cof.data());
    gemm_cpu(nc1, p.data(), cof.data(), f.data());
    if (std::getenv("PRJFC_PROBE")) {
        static bool dumped2 = false;
        if (!dumped2) {
            dumped2 = true;
            std::printf("PRJFC after-gemm f[1..8] =");
            for (int k = 1; k <= 8; ++k)
                std::printf(" %20.14f", f[k]);
            std::printf("\nPRJFC after-gemm f[%d]= %20.14f\n", nc1, f[nc1]);
            FILE* fp = std::fopen("prjfc_f_our.txt", "w");
            if (fp) {
                for (int ii = 1; ii <= nc1; ++ii)
                    for (int jj = 1; jj <= ii; ++jj)
                        std::fprintf(fp, "%.16e\n", f[jj + (ii - 1) * nc1]);
                std::fclose(fp);
            }
            fp = std::fopen("prjfc_x_our.txt", "w");
            if (fp) {
                for (int ii = 1; ii <= natm_loc; ++ii)
                    std::fprintf(fp, "%4d % .17e % .17e % .17e % .17e\n", ii,
                                 coord[0][ii], coord[1][ii], coord[2][ii],
                                 atmass[ii]);
                for (int ii = 1; ii <= 12; ++ii)
                    std::fprintf(fp, "%.17e\n", x[ii]);
                std::fprintf(fp, "ROT\n");
                for (int ii = 1; ii <= 3; ++ii)
                    std::fprintf(fp, " % .17e % .17e % .17e\n", rot[1][ii],
                                 rot[2][ii], rot[3][ii]);
                std::fclose(fp);
            }
        }
    }
}
// ---------------------------------------------------------------------------
// formd: form geometry step (NR / P-RFO / QA) (ef.F90 lines 1095-1534).
// ---------------------------------------------------------------------------
static void formd(const std::vector<double>& eigval,
                  const std::vector<double>& fx, int nvar, double dmax,
                  double ddmin, bool ts, bool& lorjk, bool rrscal, bool donr,
                  const std::vector<double>& u, std::vector<double>& d,
                  std::vector<double>& vmode) {
    const double big = 1.0e3, half = 0.5, tmtwo = 1.0e-2;
    const double eps = 1.0e-12, four = 4.0, one = 1.0, sfix = 1.0e1,
                 step_loc = 5.0e-2, ten = 1.0e1, tmsix = 1.0e-6, toll = 1.0e-8,
                 zero = 0.0;
    static int icalcn = 0;
    static double eone = 0.0, eigit = 0.0, store_ddx = 0.0;
    bool frodo1 = false, frodo2 = false, rscal;
    int i, it = 0, j, jt = 1, ncnt, newmod;
    double bl, bu, d2max, fl, fm, fu, ssmax, ssmin, sstoll, temp, xlamda,
        lamda = 0.0, lamda0 = 0.0, sstep;
    if (icalcn != numcal) {
        icalcn = numcal;
        store_ddx = 0.0;
        for (i = 1; i <= nvar; ++i) d[i] = 0.0;
    }
    skal = one;
    rscal = rrscal;
    it = 0;
    jt = 1;
    if (ts) {
        if (ef_mode != 0) {
            overlp(dmax, ddmin, newmod, nvar, lorjk, u, vmode);
            if (lorjk) return;
            if (newmod != ef_mode && iprnt >= 1)
                std::printf(" WARNING! MODE SWITCHING. WAS FOLLOWING MODE %3d "
                            "NOW FOLLOWING MODE %3d\n",
                            ef_mode, newmod);
            ef_mode = newmod;
            it = ef_mode;
        } else {
            it = 1;
        }
        eigit = eigval[it];
        if (eigit == zero) {
            if (nvar > 3 * numat - 5) {
                mopend("Too many variables. By definition, at least one force "
                       "constant is exactly zero");
                std::printf("and the lowest force constant is not negative.\n");
                std::printf("Number of variables = %5d\n", nvar);
                std::printf("Number of atoms     = %5d\n", numat);
                if (nvar == 3 * numat)
                    std::printf("(If coordinates are Cartesian, convert to "
                                "internal coordinates and re-run.)\n");
            } else {
                std::printf(" At least one force constant is exactly zero\n");
                mopend("At least one force constant is exactly zero");
            }
            return;
        } else if (iprnt >= 1) {
            std::printf("\n TS MODE IS NUMBER%3d WITH EIGENVALUE%9.1f\n"
                        "AND COMPONENTS\n",
                        it, eigit);
            for (i = 1; i <= nvar; ++i) std::printf("%9.4f", u[i + (it - 1) * nvar]);
            std::printf("\n");
        }
    }
    for (i = 1; i <= nvar; ++i) {
        if (i != it && eigval[i] != zero) {
            jt = i;
            break;
        }
    }
    eone = eigval[jt];
    ssmin = std::max(std::abs(eone) * eps, ten * eps);
    ssmax = std::max(big, std::abs(eone));
    ssmax = ssmax * big;
    sstoll = toll;
    d2max = dmax * dmax;
    frodo1 = false;
    frodo2 = false;
    lamda = zero;
    lamda0 = zero;
    if (ts) {
        if (eigit < zero && eone >= zero && donr) {
            if (iprnt >= 1)
                std::printf(" TS SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP\n");
            goto L1010;
        }
    } else if (donr) {
        if (eone >= zero) {
            if (iprnt >= 1)
                std::printf(" MIN SEARCH, CORRECT HESSIAN, TRYING PURE NR STEP\n");
            goto L1010;
        }
    }
L1000:
    if (ts) {
        lamda0 = eigval[it] + std::sqrt(eigval[it] * eigval[it] + four * fx[it] * fx[it]);
        lamda0 = lamda0 * half;
        if (iprnt >= 1)
            std::printf(" LAMDA THAT MAXIMIZES ALONG TS MODES =   %15.5f\n",
                        lamda0);
    }
    sstep = step_loc;
    if (eone <= zero) lamda = eone - sstep;
    if (eone > zero) sstep = eone;
    bl = lamda - sstep;
    bu = lamda + sstep * half;
    ncnt = 0;
    for (;;) {
        fl = zero;
        fu = zero;
        for (i = 1; i <= nvar; ++i) {
            if (i != it) {
                if ((bl - eigval[i] == zero) || (bu - eigval[i] == zero)) {
                    std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                    mopend("in FORMD");
                    return;
                }
                fl = fl + (fx[i] * fx[i]) / (bl - eigval[i]);
                fu = fu + (fx[i] * fx[i]) / (bu - eigval[i]);
            }
        }
        fl = fl - bl;
        fu = fu - bu;
        if (fl * fu < zero) break;
        bl = bl - (eone - bl);
        bu = bu + half * (eone - bu);
        if (bl <= -ssmax) {
            bl = -ssmax;
            frodo1 = true;
        }
        if (std::abs(eone - bu) <= ssmin) {
            bu = eone - ssmin;
            frodo2 = true;
        }
        if (frodo1 && frodo2) {
            std::printf("NUMERICAL PROBLEMS IN BRACKETING LAMDA %g %g %g %g %g\n",
                        eone, bl, bu, fl, fu);
            std::printf(" GOING FOR FIXED STEP SIZE....\n");
            goto L1020;
        }
        ++ncnt;
        if (ncnt > 1000) {
            std::printf("TOO MANY ITERATIONS IN LAMDA BISECT %g %g %g %g\n",
                        bl, bu, fl, fu);
            mopend("TOO MANY ITERATIONS IN LAMDA BISECT IN EF");
            return;
        }
    }
    ncnt = 0;
    xlamda = zero;
    for (;;) {
        fl = zero;
        fu = zero;
        fm = zero;
        lamda = half * (bl + bu);
        for (i = 1; i <= nvar; ++i) {
            if (i != it) {
                if ((bl - eigval[i] == zero) || (bu - eigval[i] == zero) ||
                    (lamda - eigval[i] == zero)) {
                    std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                    mopend("in FORMD");
                    return;
                }
                fl = fl + (fx[i] * fx[i]) / (bl - eigval[i]);
                fu = fu + (fx[i] * fx[i]) / (bu - eigval[i]);
                fm = fm + (fx[i] * fx[i]) / (lamda - eigval[i]);
            }
        }
        fl = fl - bl;
        fu = fu - bu;
        fm = fm - lamda;
        if (std::abs(xlamda - lamda) < sstoll) break;
        ++ncnt;
        if (ncnt > 1000) {
            std::printf("TOO MANY ITERATIONS IN LAMDA BISECT %g %g %g %g %g\n",
                        bl, bu, lamda, fl, fu);
            mopend("TOO MANY ITERATIONS IN LAMDA BISECT IN EF");
            return;
        }
        xlamda = lamda;
        if (fm * fu < zero) bl = lamda;
        if (fm * fl < zero) bu = lamda;
    }
L1010:
    if (iprnt >= 1)
        std::printf(" LAMDA THAT MINIMIZES ALONG ALL MODES =  %17.9f\n", lamda);
    for (i = 1; i <= nvar; ++i) d[i] = zero;
    for (i = 1; i <= nvar; ++i) {
        if (lamda == zero && std::abs(eigval[i]) < tmtwo)
            temp = zero;
        else
            temp = fx[i] / (lamda - eigval[i]);
        if (i == it) {
            if (std::abs(lamda0 - eigval[it]) < 1.0e-9) {
                std::printf(" TS FAILED TO LOCATE TRANSITION STATE\n");
                std::printf("\n CURRENT VALUE OF geo\n");
                geout(iw);
                mopend("TS FAILED TO LOCATE TRANSITION STATE");
                return;
            } else {
                temp = fx[it] / (lamda0 - eigval[it]);
            }
        }
        if (iprnt >= 5) std::printf("FORMD, DELTA STEP %d %16.11f\n", i, temp);
        for (j = 1; j <= nvar; ++j) d[j] = d[j] + temp * u[j + (i - 1) * nvar];
    }
    ddx = std::sqrt(ddot(nvar, &d[1], 1, &d[1], 1));
    if ((lamda != zero && std::abs(lamda0 + lamda) <= 1.0e-20) ||
        (lamda != zero && std::abs(lamda0 + lamda) > 1.0e-20 &&
         store_ddx > 1.0e-9)) {
        if (ddx > 2 * store_ddx) {
            for (i = 1; i <= nvar; ++i) d[i] = d[i] * store_ddx / ddx;
            ddx = 2 * store_ddx;
        }
    }
    if (lamda == zero && lamda0 == zero && iprnt >= 1)
        std::printf(" PURE NR-STEP HAS LENGTH%10.5f\n", ddx);
    if (lamda != zero && std::abs(lamda0 + lamda) > 1.0e-20 && iprnt >= 1)
        std::printf(" P-RFO-STEP   HAS LENGTH%10.5f\n", ddx);
    if (lamda != zero && std::abs(lamda0 + lamda) <= 1.0e-20 && iprnt >= 1)
        std::printf(" QA/TRIM-STEP HAS LENGTH%10.5f\n", ddx);
    store_ddx = ddx;
    if (ddx < (dmax + tmsix)) {
        xlamd = lamda;
        xlamd0 = lamda0;
        return;
    }
    if (lamda == zero && lamda0 == zero) goto L1000;
    if (rscal) goto L1040;
L1020:
    lamda = zero;
    frodo1 = false;
    frodo2 = false;
    sstep = step_loc;
    if (eone <= zero) lamda = eone - sstep;
    if (ts && -eigit < eone) lamda = -eigit - sstep;
    if (eone > zero) sstep = eone;
    bl = lamda - sstep;
    bu = lamda + sstep * half;
    for (;;) {
        fl = zero;
        fu = zero;
        for (i = 1; i <= nvar; ++i) {
            if (i != it) {
                if ((bl - eigval[i] == zero) || (bu - eigval[i] == zero)) {
                    std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                    mopend("in FORMD");
                    return;
                }
                double t1 = fx[i] / (bl - eigval[i]);
                double t2 = fx[i] / (bu - eigval[i]);
                fl = fl + t1 * t1;
                fu = fu + t2 * t2;
            }
        }
        if (ts) {
            if ((bl + eigval[it] == zero) || (bu + eigval[it] == zero)) {
                std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                mopend("in FORMD");
                return;
            }
            double t1 = fx[it] / (bl + eigval[it]);
            double t2 = fx[it] / (bu + eigval[it]);
            fl = fl + t1 * t1;
            fu = fu + t2 * t2;
        }
        fl = fl - d2max;
        fu = fu - d2max;
        if (fl * fu < zero) goto L1030;
        bl = bl - (eone - bl);
        bu = bu + half * (eone - bu);
        if (bl <= -ssmax) {
            bl = -ssmax;
            frodo1 = true;
        }
        if (std::abs(eone - bu) <= ssmin) {
            bu = eone - ssmin;
            frodo2 = true;
        }
        if (frodo1 && frodo2) break;
    }
    std::printf("NUMERICAL PROBLEMS IN BRACKETING LAMDA %g %g %g %g %g\n",
                eone, bl, bu, fl, fu);
    std::printf(" GOING FOR FIXED LEVEL SHIFTED NR STEP...\n");
    lamda = eone - sfix;
    lamda0 = eigit + sfix;
    rscal = true;
    goto L1010;
L1030:
    ncnt = 0;
    xlamda = zero;
    for (;;) {
        fl = zero;
        fu = zero;
        fm = zero;
        lamda = half * (bl + bu);
        for (i = 1; i <= nvar; ++i) {
            if (i != it) {
                if ((bl - eigval[i] == zero) || (bu - eigval[i] == zero) ||
                    (lamda - eigval[i] == zero)) {
                    std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                    mopend("in FORMD");
                    return;
                }
                double t1 = fx[i] / (bl - eigval[i]);
                double t2 = fx[i] / (bu - eigval[i]);
                double t3 = fx[i] / (lamda - eigval[i]);
                fl = fl + t1 * t1;
                fu = fu + t2 * t2;
                fm = fm + t3 * t3;
            }
        }
        if (ts) {
            if ((bl + eigval[it] == zero) || (bu + eigval[it] == zero) ||
                (lamda + eigval[it] == zero)) {
                std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                mopend("in FORMD");
                return;
            }
            double t1 = fx[it] / (bl + eigval[it]);
            double t2 = fx[it] / (bu + eigval[it]);
            double t3 = fx[it] / (lamda + eigval[it]);
            fl = fl + t1 * t1;
            fu = fu + t2 * t2;
            fm = fm + t3 * t3;
        }
        fl = fl - d2max;
        fu = fu - d2max;
        fm = fm - d2max;
        if (std::abs(xlamda - lamda) < sstoll) break;
        ++ncnt;
        if (ncnt > 1000) {
            std::printf("TOO MANY ITERATIONS IN LAMDA BISECT %g %g %g %g %g\n",
                        bl, bu, lamda, fl, fu);
            mopend("TOO MANY ITERATIONS IN LAMDA BISECT IN EF");
            return;
        }
        xlamda = lamda;
        if (fm * fu < zero) bl = lamda;
        if (fm * fl < zero) bu = lamda;
    }
    lamda0 = -lamda;
    rscal = true;
    goto L1010;
L1040:
    if (ddx == zero) {
        std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
        mopend("in FORMD");
        return;
    }
    skal = dmax / ddx;
    for (i = 1; i <= nvar; ++i) d[i] = d[i] * skal;
    ddx = std::sqrt(ddot(nvar, &d[1], 1, &d[1], 1));
    if (iprnt >= 1)
        std::printf(" CALCULATED STEP SIZE TOO LARGE, SCALED WITH%9.5f\n", skal);
    xlamd = lamda;
    xlamd0 = lamda0;
}
// ---------------------------------------------------------------------------
// ef: main EF optimizer driver (ef.F90 lines 1-740).
// ---------------------------------------------------------------------------
void ef(std::vector<double>& xparam, double& funct) {
    static int icalcn = 0, igthes = 0, ihess, iloop, ireclc, iupd = 0, ldump,
               mxstep, ntime;
    bool lorjk, lrjk, lts, old, let, lupd, scf1, newhes, saddle, rrscal,
        first_time, donr;
    int i, ij, imode, instep, itry1, ittest, j, k, l, neg, nflush, ih;
    double ss, absmin, deact, depre, dtmp, olde, ratio, tprt, tstep, tt0, xtmp,
        rmx, best_funct, best_gnorm, cosine, cos_const;
    double time1, time2, tol2, dmax, ddmin, ddmax, t_hess1, t_hess2;
    double odmax = 0.0, oolde = 0.0, odd = 0.0;
    char txt;
    std::vector<int> ipow(10, 0);
    std::vector<int> best_nc;
    std::vector<double> pmat, bmat, u, eigval, tvec, svec, fx, oldfx,
        oldeig, ooldf, oldf, d, vmode, best_xparam, best_grad;
    ih = (nvar * (nvar + 1)) / 2;
    if (nvar > 7996) {
        std::printf("Insufficient memory to run EF\n");
        std::printf("The L-BFGS geometry optimizer uses less memory\n"
                    "To use the L-BFGS geometry optimizer, add LBFGS to the "
                    "keyword line\n");
        mopend("Insufficient memory to run EF");
        return;
    }
    pmat.assign(nvar * nvar + 1, 0.0);
    bmat.assign(nvar * nvar + 1, 0.0);
    u.assign(nvar * nvar + 1, 0.0);
    eigval.assign(nvar + 1, 0.0);
    tvec.assign(nvar + 1, 0.0);
    svec.assign(nvar + 1, 0.0);
    fx.assign(nvar + 1, 0.0);
    oldfx.assign(nvar + 1, 0.0);
    oldeig.assign(nvar + 1, 0.0);
    ooldf.assign(nvar + 1, 0.0);
    oldf.assign(nvar + 1, 0.0);
    d.assign(nvar + 1, 0.0);
    vmode.assign(nvar + 1, 0.0);
    best_xparam.assign(nvar + 1, 0.0);
    best_grad.assign(nvar + 1, 0.0);
    best_nc.assign(natoms + 1, 0);
    best_funct = 1.0e10;
    best_gnorm = 1.0e10;
    cos_const = 1.0;
    ihess = 0;
    nstep = 0;
    for (i = 1; i <= 9; ++i) ipow[i] = 0;
    int isz = (int)hesinv.size() > 0 ? (int)hesinv.size() - 1 : 0;
    if (icalcn != numcal || nvar * nvar != isz) {
        newhes = true;
        old = (keywrd.find(" OLD_HESS") != std::string::npos);
        if (old) {
            old = (int)hesinv.size() > 0;
            if (!old) {
                line = " OLD_HESS requested, but old Hessian is missing";
                mopend(line.c_str());
            }
            if (old) {
                old = (nvar * nvar == isz);
                if (!old) {
                    line = " OLD_HESS requested, but old Hessian is of a different size";
                    mopend(line.c_str());
                }
            }
            if (!old) return;
        }
        if (!old) {
            hesinv.assign(nvar * nvar + 1, 0.0);
        }
        efstr(xparam, funct, ihess, ntime, iloop, igthes, mxstep, ireclc, iupd,
              dmax, ddmax, ddmin, tol2, time1, time2, nvar, scf1, lupd, ldump,
              rrscal, donr, hesinv, bmat, pmat, grad, oldf, d, vmode);
        if (moperr) goto L1100;
        if (old) iloop = -1;
        let = (keywrd.find(" LET") != std::string::npos);
        saddle = (keywrd.find(" SADDLE") != std::string::npos);
        first_time = true;
    } else {
        first_time = false;
    }
    last = 1;
    for (i = 1; i <= nvar; ++i) grad[i] = 0.0;
    time1 = second(1);
    compfg(xparam, true, funct, true, grad, true);
    if (moperr) goto L1100;
    olde = funct;
    nflush = 1;
    absmin = 1.0e9;
    lts = false;
    if (negreq == 1) lts = true;
    lorjk = false;
    if (scf1) {
        gnorm = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        iflepo = 1;
        goto L1100;
    }
    rmx = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
    if (rmx < tol2) {
        iflepo = 2;
        last = 1;
        icalcn = numcal;
        std::printf("\n GRADIENT =%9.5f  IS LESS THAN CUTOFF =%9.5f\n\n", rmx,
                    tol2);
        goto L1030;
    }
    last = 0;
    t_hess1 = second(1);
    if (newhes && iloop > 0) {
        gethes(xparam, igthes, iloop, hesinv, pmat, bmat, grad, geo, loc, oldf,
               d, vmode, funct);
        if (moperr) goto L1100;
        newhes = false;
    }
    t_hess2 = second(1);
    icalcn = numcal;
    iflepo = 0;
    funct = olde;
    for (;;) {
        for (i = 1; i <= nvar; ++i) {
            oldfx[i] = fx[i];
            ooldf[i] = oldf[i];
            oldeig[i] = eigval[i];
            for (j = 1; j <= nvar; ++j) {
                bmat[i + (j - 1) * nvar] = hesinv[i + (j - 1) * nvar];
                pmat[i + (j - 1) * nvar] = u[i + (j - 1) * nvar];
            }
        }
        if (ihess >= ireclc && iflepo != 15) {
            iloop = 1;
            ihess = 0;
            if (igthes != 3) igthes = 1;
            gethes(xparam, igthes, iloop, hesinv, pmat, bmat, grad, geo, loc,
                   oldf, d, vmode, funct);
            if (moperr) goto L1100;
        }
        if (ihess > 0) {
            updhes(svec, tvec, grad, nvar, iupd, hesinv, oldf, d);
        }
        if (iprnt >= 2) geout(iw);
        fprintf(stderr, "[EFGRAD] nstep=%d grad=%+.12f %+.12f %+.12f\n", nstep, grad[1], grad[2], grad[3]);
        if (iprnt >= 2) {
            std::printf(" XPARAM \n");
            for (i = 1; i <= nvar; ++i)
                std::printf(" %2d%3d%10.4f\n", loc[1][i], loc[2][i], xparam[i]);
            std::printf(" grad\n");
            for (i = 1; i <= nvar; ++i) std::printf("%11.5f", grad[i]);
            std::printf("\n");
        }
        gnorm = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        time2 = second(2);
        tstep = time2 - time1;
        if (tstep < 0.0) tstep = 0.0;
        tleft = tleft - tstep;
        if (tleft < 0.0) tleft = -0.1;
        time1 = time2;
        if (tleft < (tstep - (t_hess2 - t_hess1)) * 2.0) goto L1030;
        if ((!saddle && tstep > 0.5) || first_time) {
            prttim(tleft, tprt, txt);
            if (ldump == 0) {
                if (id == 3) {
                    write_cell(iw);
                    write_cell(iw0);
                }
                std::printf(" CYCLE:%6d TIME:%8.3f TIME LEFT:%6.2f%c  "
                            "GRAD.:%10.3f HEAT:%14.7g\n",
                            nstep + 1, std::min(tstep, 9999.99), tprt, txt,
                            std::min(gnorm, 999999.999), funct);
                to_screen("To_file: Geometry optimizing");
            } else {
                if (id == 3) {
                    write_cell(iw);
                    write_cell(iw0);
                }
                std::printf(" RESTART FILE WRITTEN,      TIME LEFT:%6.2f%c  "
                            "GRAD.:%10.3f HEAT:%14.7g\n",
                            tprt, txt, std::min(gnorm, 999999.999), funct);
                to_screen("To_file: Geometry optimizing");
            }
        }
        if ((ef_mode != 0 && best_gnorm > gnorm) ||
            (ef_mode == 0 && best_funct > funct)) {
            best_gnorm = gnorm;
            best_funct = funct;
            for (i = 1; i <= nvar; ++i) {
                best_xparam[i] = xparam[i];
                best_grad[i] = grad[i];
            }
            for (i = 1; i <= natoms; ++i) best_nc[i] = nc[i];
        }
        // label 1000 in Fortran
        ihess = ihess + 1;
        nstep = nstep + 1;
        rmx = std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1));
        if (rmx < tol2) {
            std::printf("\n GRADIENT =%9.5f  IS LESS THAN CUTOFF =%9.5f\n\n",
                        rmx, tol2);
            goto L1020;
        } else {
            if (ef_mode == 0 && nstep > 5) {
                if (absmin - funct < 1.0e-5) {
                    if (funct - absmin > 1.0) itry1 = 0;
                    if (itry1 > 900 || (rmx < 0.1 && itry1 > 9)) {
                        std::printf(
                            "\n HEAT OF FORMATION IS ESSENTIALLY STATIONARY\n");
                        goto L1020;
                    }
                    itry1 = itry1 + 1;
                } else {
                    itry1 = 0;
                    absmin = funct;
                }
            }
            olde = funct;
            if (nstep > 5) {
                cosine = ddot(nvar, &grad[1], 1, &oldf[1], 1) /
                         std::sqrt(ddot(nvar, &grad[1], 1, &grad[1], 1) *
                                   ddot(nvar, &oldf[1], 1, &oldf[1], 1));
                if (cosine > 0.9)
                    cos_const = cos_const * 1.5;
                else
                    cos_const = 1.0;
            }
            for (i = 1; i <= nvar; ++i) oldf[i] = grad[i];
            if (nvar == 3 * numat && numat > 2) {
                // Static buffers: prjfc results are consumed only inside prjfc's
                // call (projection side-effect), so reuse avoids 5 heap
                // allocations per EF iteration (78^2 doubles + 4x coord_l).
                static std::vector<double> p_l, x_l, rm_l, dx_l;
                static std::vector<std::vector<double>> coord_l;
                if ((int)p_l.size() < nvar * nvar + 1) {
                    p_l.assign(nvar * nvar + 1, 0.0);
                    x_l.assign(nvar + 1, 0.0);
                    rm_l.assign(nvar + 1, 0.0);
                    dx_l.assign(nvar + 1, 0.0);
                    coord_l.assign(4, std::vector<double>(numat + 1, 0.0));
                }
                prjfc(hesinv, xparam, nvar, u, p_l, atmass, x_l, rm_l, dx_l,
                      coord_l);
                if (moperr) goto L1100;
            }
            {
                std::vector<double> rsp_in(nvar * nvar, 0.0),
                    rsp_root(nvar, 0.0), rsp_vec(nvar * nvar, 0.0);
                ij = 0;
                for (i = 1; i <= nvar; ++i)
                    for (j = 1; j <= i; ++j) {
                        rsp_in[ij] = hesinv[j + (i - 1) * nvar];
                        ++ij;
                    }
                if (std::getenv("HESSPROBE")) {
                    std::printf(" HESSPROBE cyc%6d\n", nstep + 1);
                    std::printf("  hc ");
                    int ihp = nvar * (nvar + 1) / 2;
                    for (int k = 0; k < ihp; ++k) std::printf("%20.14f", rsp_in[k]);
                    std::printf("\n");
                }
                // ULP perturbation injection (research probe): MOPAC_ULPPERT=idx:dir
                // idx: 1-based linear index into packed hessc (1..nvar*(nvar+1)/2);
                // dir: +1/-1 (nextafter). Applied once on the first EF cycle only.
                {
                    static bool ulp_done = false;
                    const char* upp = std::getenv("MOPAC_ULPPERT");
                    if (upp && !ulp_done) {
                        int uidx = 0, udir = 1;
                        if (std::sscanf(upp, "%d:%d", &uidx, &udir) == 2 &&
                            uidx >= 1 && uidx <= nvar * (nvar + 1) / 2) {
                            double base = rsp_in[uidx - 1];
                            rsp_in[uidx - 1] = std::nextafter(
                                base, (udir > 0) ? std::numeric_limits<double>::infinity()
                                                 : -std::numeric_limits<double>::infinity());
                            std::printf(" ULPPERT idx=%d dir=%d base=%.17e new=%.17e\n",
                                        uidx, udir, base, rsp_in[uidx - 1]);
                            ulp_done = true;
                        }
                    }
                }
            rsp(rsp_in.data(), nvar, rsp_root.data(), rsp_vec.data());
                for (i = 1; i <= nvar; ++i) {
                    eigval[i] = rsp_root[i - 1];
                    for (j = 1; j <= nvar; ++j)
                        u[j + (i - 1) * nvar] =
                            rsp_vec[(i - 1) * nvar + (j - 1)];
                }
            }
            for (i = 1; i <= nvar; ++i) {
                if (std::abs(eigval[i]) < 1.0e-6) eigval[i] = 0.0;
            }
            ij = nvar * nvar;
            for (i = nvar; i >= 1; --i)
                for (j = nvar; j >= 1; --j) {
                    u[j + (i - 1) * nvar] = u[ij];
                    --ij;
                }
            if (iprnt >= 3) prthes(eigval, nvar, hesinv, u);
            if (mxstep == 0) {
                nstep = 0;
                goto L1030;
            }
            neg = 0;
            for (i = 1; i <= nvar; ++i) {
                if (eigval[i] < 0.0) neg = neg + 1;
            }
            if (iprnt >= 1) {
                std::printf("\n HESSIAN HAS%3d NEGATIVE EIGENVALUE(S)\n", neg);
                for (i = 1; i <= neg; ++i) std::printf("%7.1f", eigval[i]);
                std::printf("\n");
            }
            // fx = row_i(u) . grad  (official Fortran semantics; verified
            // bit-identical to the official exe on the 10 small molecules and
            // required for EF convergence — column-form/skip-degenerate
            // variants diverge, see paper-discussion notes 2026-10-06).
            // Optional DSPG variant (MOPAC_DSPG=1, paper path C): project grad
            // onto the non-degenerate subspace (grad' = grad - P_d*grad, P_d
            // the |eigval|<1e-6-subspace projector) before computing fx.
            // C12H26: 267 steps / -66.77024 vs baseline 227 / -66.52132 —
            // under test on more molecules for systematic improvement.
            // ---- GAB: gradient-aligned degenerate basis (research probe) ----
            // Rotate the degenerate columns of u so that the first one aligns
            // with the gradient projection onto the degenerate subspace.
            // This concentrates fx into one degenerate component, giving EF
            // clean per-mode decisions (hypothesis: official dsyevd basis
            // performs well for the same reason).  MOPAC_GAB=1 gated, default
            // off (official path stays bit-identical).
            if (std::getenv("MOPAC_GAB")) {
                std::vector<int> D;
                for (i = 1; i <= nvar; ++i)
                    if (std::abs(eigval[i]) < 1.0e-6) D.push_back(i);
                int m = (int)D.size();
                if (m >= 2 && m < nvar) {
                    std::vector<double> gd(m, 0.0);
                    for (int a = 0; a < m; ++a) {
                        int col = D[a];
                        for (int r = 1; r <= nvar; ++r)
                            gd[a] += u[r + (col - 1) * nvar] * grad[r];
                    }
                    double ng = 0.0;
                    for (int a = 0; a < m; ++a) ng += gd[a] * gd[a];
                    ng = std::sqrt(ng);
                    if (ng > 0.0) {
                        std::vector<std::vector<double>> Q(m,
                                                           std::vector<double>(m, 0.0));
                        for (int a = 0; a < m; ++a) Q[0][a] = gd[a] / ng;
                        for (int b = 1; b < m; ++b) {
                            for (int a = 0; a < m; ++a)
                                Q[b][a] = (a == b - 1) ? 1.0 : 0.0;
                            for (int c = 0; c < b; ++c) {
                                double dot = 0.0;
                                for (int a = 0; a < m; ++a)
                                    dot += Q[b][a] * Q[c][a];
                                for (int a = 0; a < m; ++a)
                                    Q[b][a] -= dot * Q[c][a];
                            }
                            double nn = 0.0;
                            for (int a = 0; a < m; ++a) nn += Q[b][a] * Q[b][a];
                            if (nn > 0.0) {
                                nn = std::sqrt(nn);
                                for (int a = 0; a < m; ++a) Q[b][a] /= nn;
                            } else {
                                for (int a = 0; a < m; ++a)
                                    Q[b][a] = (a == b) ? 1.0 : 0.0;
                            }
                        }
                        std::vector<double> tmp(nvar * m, 0.0);
                        for (int b = 0; b < m; ++b)
                            for (int r = 1; r <= nvar; ++r) {
                                double s = 0.0;
                                for (int a = 0; a < m; ++a)
                                    s += Q[b][a] * u[r + (D[a] - 1) * nvar];
                                tmp[b * nvar + (r - 1)] = s;
                            }
                        for (int b = 0; b < m; ++b)
                            for (int r = 1; r <= nvar; ++r)
                                u[r + (D[b] - 1) * nvar] =
                                    tmp[b * nvar + (r - 1)];
                        std::printf(" GAB m=%d |gd|=%.6e fx-concentrated\n",
                                    m, ng);
                    }
                }
            }
            if (std::getenv("MOPAC_DSPG")) {
                std::vector<double> gradp(nvar + 1, 0.0);
                for (i = 1; i <= nvar; ++i) gradp[i] = grad[i];
                for (k = 1; k <= nvar; ++k) {
                    if (std::abs(eigval[k]) == 0.0) {
                        double c =
                            ddot(nvar, u.data() + 1 + (k - 1) * nvar, 1,
                                 &grad[1], 1);
                        for (int r = 1; r <= nvar; ++r)
                            gradp[r] -= c * u[r + (k - 1) * nvar];
                    }
                }
                for (i = 1; i <= nvar; ++i) {
                    fx[i] = ddot(nvar, u.data() + 1 + (i - 1) * nvar, 1,
                                 &gradp[1], 1);
                    if (std::abs(eigval[i]) == 0.0) fx[i] = 0.0;
                }
            } else {
                for (i = 1; i <= nvar; ++i) {
                    fx[i] = ddot(nvar, u.data() + 1 + (i - 1) * nvar, 1,
                                 &grad[1], 1);
                    if (std::abs(eigval[i]) == 0.0) fx[i] = 0.0;
                }
            }
            if (ihess == 1 && std::getenv("EFPROBE")) {
                std::printf("\n  fc ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", fx[i]);
                std::printf("\n  ec ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", eigval[i]);
                std::printf("\n  u8c ");
                for (j = 1; j <= 5; ++j) std::printf("%13.7f", u[j + 7 * nvar]);
                std::printf("\n  u7v ");
                for (j = 1; j <= 5; ++j)
                    std::printf("%13.7f", u[7 + (j - 1) * nvar]);
                std::printf("\n");
            }
            // EFPROBE: diagnostics only (env EFPROBE=1), mirrors official
            // ef.F90 probe block — prints xparam/grad/eigval/fx/d per cycle.
            if (std::getenv("EFPROBE")) {
                std::printf(" EFPROBE cyc%6d\n", nstep + 1);
                std::printf("  xp ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", xparam[i]);
                std::printf("\n  gr ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", grad[i]);
                std::printf("\n  ev ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", eigval[i]);
                std::printf("\n  fx ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", fx[i]);
                std::printf("\n  d  ");
                for (i = 1; i <= nvar; ++i) std::printf("%13.7f", d[i]);
                std::printf("\n  u7 ");
                for (j = 1; j <= 10; ++j) std::printf("%13.7f", u[j + 6 * nvar]);
                std::printf("\n  u8 ");
                for (j = 1; j <= 10; ++j) std::printf("%13.7f", u[j + 7 * nvar]);
                std::printf("\n");
            }
            for (;;) {
                formd(eigval, fx, nvar, dmax, ddmin, lts, lorjk, rrscal, donr,
                      u, d, vmode);
                if (std::getenv("EFPROBE")) {
                    std::printf("\n  dn ");
                    for (i = 1; i <= nvar; ++i) std::printf("%13.7f", d[i]);
                    std::printf("\n");
                }
                if (moperr) goto L1100;
                if (lorjk) {
                    if (iprnt >= 1) std::printf("      NOW UNDOING PREVIOUS STEP\n");
                    dmax = odmax;
                    ddx = odd;
                    olde = oolde;
                    for (i = 1; i <= nvar; ++i) {
                        fx[i] = oldfx[i];
                        oldf[i] = ooldf[i];
                        eigval[i] = oldeig[i];
                        for (j = 1; j <= nvar; ++j) {
                            hesinv[i + (j - 1) * nvar] = bmat[i + (j - 1) * nvar];
                            u[i + (j - 1) * nvar] = pmat[i + (j - 1) * nvar];
                        }
                    }
                    for (i = 1; i <= nvar; ++i) {
                        xparam[i] = xparam[i] - d[i];
                        k = loc[1][i];
                        l = loc[2][i];
                        geo[l][k] = xparam[i];
                    }
                    if (ndep != 0) symtry();
                    dmax = std::min(dmax, ddx) / 2.0;
                    odmax = dmax;
                    odd = ddx;
                    nstep = nstep - 1;
                    if (dmax < ddmin) {
                        goto L1010;
                    } else if (iprnt >= 1) {
                        std::printf("      FINISH UNDOING, NOW GOING FOR NEW STEP\n");
                    }
                } else {
                    for (i = 1; i <= nvar; ++i) d[i] = d[i] * cos_const;
                    fprintf(stderr, "[EFSTEP] nstep=%d cos=%f d=%+.12f %+.12f %+.12f\n", nstep, cos_const, d[1], d[2], d[3]);
                    for (i = 1; i <= nvar; ++i) {
                        xparam[i] = xparam[i] + d[i];
                        k = loc[1][i];
                        l = loc[2][i];
                        geo[l][k] = xparam[i];
                    }
                    if (ndep != 0) symtry();
                    depre = 0.0;
                    imode = 1;
                    if (ef_mode != 0) imode = ef_mode;
                    for (i = 1; i <= nvar; ++i) {
                        if (lts && i == imode)
                            xtmp = xlamd0;
                        else
                            xtmp = xlamd;
                        if (std::abs(xtmp - eigval[i]) >= 1.0e-2) {
                            ss = skal * fx[i] / (xtmp - eigval[i]);
                            depre = depre + ss * (fx[i] + 0.5 * ss * eigval[i]);
                        }
                    }
                    for (i = 1; i <= nvar; ++i) grad[i] = 0.0;
                    compfg(xparam, true, funct, true, grad, true);
                    if (moperr) goto L1100;
                    deact = funct - olde;
                    if (depre == 0.0) {
                        std::printf(" CALCULATION IS TERMINATED TO AVOID ZERO DIVIDE\n");
                        mopend("in EF");
                        goto L1100;
                    }
                    ratio = deact / depre;
                    if (iprnt >= 1) {
                        std::printf("       HoF         ACTUAL,  PREDICTED "
                                    "ENERGY CHANGE, RATIO\n %14.7f%20.7f%13.7f%12.4f\n",
                                    funct, deact, depre, ratio);
                    }
                    lrjk = false;
                    if ((ratio < rmin || ratio > rmax) &&
                        (std::abs(depre) > 1.0e-2 || std::abs(deact) > 1.0e-2)) {
                        dtmp = std::min(dmax, ddx) / 2.0;
                        if (dtmp <= ddmin) dtmp = ddmin;
                        if (iprnt >= 1)
                            std::printf("UNACCEPTABLE RATIO, REJECTING STEP, "
                                        "REDUCING DMAX TO%7.4f\n",
                                        dtmp);
                        lrjk = true;
                    }
                    if (lrjk && std::abs(dmax - ddmin) < 1.0e-20) {
                        if (iprnt >= 1)
                            std::printf("NEW TRUST RADIUS WOULD BE BELOW DDMIN "
                                        "ACCEPTING STEP ANYWAY\n");
                        lrjk = false;
                    }
                    if (let) lrjk = false;
                    if (!lrjk) break;
                    for (i = 1; i <= nvar; ++i) {
                        xparam[i] = xparam[i] - d[i];
                        k = loc[1][i];
                        l = loc[2][i];
                        geo[l][k] = xparam[i];
                    }
                    if (ndep != 0) symtry();
                    dmax = std::min(dmax, ddx) / 2.0;
                    if (dmax < ddmin) goto L1010;
                }
            }
            if (iprnt >= 1) std::printf(" STEPSIZE USED IS%9.5f\n", ddx);
            if (iprnt >= 2) {
            std::printf(" CALCULATED STEP\n");
                for (i = 1; i <= nvar; ++i) std::printf("%10.6f", d[i]);
                std::printf("\n");
            }
            odmax = dmax;
            odd = ddx;
            oolde = olde;
            if (lupd && ((rmx > gmin) ||
                         (std::abs(depre) > 1.0e-2 || std::abs(deact) > 1.0e-2))) {
                if (lts && (ratio <= 0.1 || ratio >= 3.0))
                    dmax = std::min(dmax, ddx) / 2.0;
                if (lts && ratio >= 0.75 && ratio <= (4.0 / 3.0) &&
                    ddx > (dmax - 1.0e-6))
                    dmax = dmax * std::sqrt(2.0);
                if (!lts && ratio >= 0.5 && ddx > (dmax - 1.0e-6))
                    dmax = dmax * std::sqrt(2.0);
                if (std::abs(ratio - 1.0) < 0.1) dmax = dmax * std::sqrt(2.0);
                dmax = std::max(dmax, ddmin);
                dmax = std::min(dmax, ddmax);
            }
            if (lupd && rmx < gmin &&
                (std::abs(depre) < 1.0e-2 && std::abs(deact) < 1.0e-2)) {
                dmax = std::max(dmax, 0.1);
            }
            if (iprnt >= 1) std::printf(" NEW TRUST RADIUS = %8.5f\n", dmax);
        }
    L1010:
        if (dmax < ddmin) break;
        if (nstep >= mxstep) goto L1030;
        ittest = (int)((time2 - time0) / tdump);
        if (ittest > ntime) {
            ldump = 1;
            ntime = std::max(ittest, (ntime + 1));
            ipow[1] = ihess;
            ipow[2] = nstep;
            ipow[9] = 2;
            tt0 = second(1) - time0;
            instep = -nstep;
            efsav(tt0, hesinv, funct, grad, xparam, pmat, instep, bmat, ipow,
                  oldf, d, vmode);
            if (moperr) goto L1100;
        } else {
            ldump = 0;
        }
    }
    std::printf("\n TRUST RADIUS NOW LESS THAN %8.5f OPTIMIZATION TERMINATING\n"
                " THE GEOMETRY MAY NOT BE COMPLETELY OPTIMIZED\n"
                " (TO CONTINUE, ADD 'LET DDMIN=0.0' TO THE KEYWORD LINE)\n",
                ddmin);
    if ((ef_mode != 0 && best_gnorm < gnorm) ||
        (ef_mode == 0 && best_funct < funct)) {
        funct = best_funct;
        gnorm = best_gnorm;
        for (i = 1; i <= nvar; ++i) {
            xparam[i] = best_xparam[i];
            grad[i] = best_grad[i];
        }
        for (i = 1; i <= natoms; ++i) nc[i] = best_nc[i];
    }
L1020:
    iflepo = 15;
    last = 1;
    tt0 = second(1) - time0;
    if (!moperr) {
        compfg(xparam, true, funct, true, grad, false);
    }
    goto L1100;
L1030:
    if (tleft < tstep * 2.0) mopend("NOT ENOUGH TIME FOR ANOTHER CYCLE");
    if (nstep >= mxstep) mopend("EXCESS NUMBER OF OPTIMIZATION CYCLES");
    if (tleft < tstep * 2.0 || nstep >= mxstep) {
        ipow[1] = ihess;
        ipow[9] = 1;
        ipow[2] = nstep;
        tt0 = second(1) - time0;
        instep = -nstep;
        efsav(tt0, hesinv, funct, grad, xparam, pmat, instep, bmat, ipow, oldf,
              d, vmode);
        iflepo = -1;
    }
L1100:
    return;
}
