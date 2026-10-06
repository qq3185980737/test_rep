// pulay.cpp — Pulay DIIS convergence.
#include "pulay.h"
#include <cmath>
#include <algorithm>
#include "molkst_C.h"
using molkst_C::numcal;
using molkst_C::keywrd;
using molkst_C::mpack;
extern void mamult(const double*, const double*, double*, int, double);
extern void osinv(double*, int, double&);
static double ddot(int n, const double* a, const double* b) {
    double s=0; for(int i=0;i<n;++i) s+=a[i]*b[i]; return s;
}

void pulay(std::vector<double>& f, std::vector<double>& p, int n,
           std::vector<double>& fppf, std::vector<double>& fock,
           std::vector<double>& diag, int& nff, int& nfock, int mdf,
           bool& reset, double& rms) {
    pulay(f.data(), p.data(), n, fppf.data(), fock.data(), diag.data(), nff, nfock, mdf, reset, rms);
}
void pulay(double* f, double* p, int n, double* fppf, double* fock,
           double* emat, int& lfock, int& nfock, int msize, bool& start,
           double& pl) {
    static int icalcn=0, maxlim=6;
    static bool debug=false;
    static int linear=0, mfock=0;
    if (icalcn != numcal) {
        icalcn = numcal;
        maxlim = 6;
        debug = keywrd.find("DEBUGPULAY") != std::string::npos;
    }
    if (start) {
        linear = n*(n+1)/2;
        mfock = msize/linear;
        mfock = std::min(maxlim, mfock);
        nfock = 1; lfock = 1; start = false;
    } else {
        if (nfock < mfock) ++nfock;
        if (lfock != mfock) ++lfock; else lfock = 1;
    }
    int lbase = (lfock-1)*linear;
    for (int i=0;i<linear;++i) fock[lfock-1 + i*mfock] = f[i];
    mamult(p, f, fppf+lbase, n, 0.0);
    mamult(f, p, fppf+lbase, n, -1.0);
    int nfock1 = nfock + 1;
    for (int i=1;i<=nfock;++i) {
        emat[(nfock1-1)*20 + i-1] = -1.0;
        emat[(i-1)*20 + nfock1-1] = -1.0;
        double v = ddot(linear, fppf+(i-1)*linear, fppf+lbase);
        emat[(lfock-1)*20 + i-1] = v;
        emat[(i-1)*20 + lfock-1] = v;
    }
    pl = emat[(lfock-1)*20 + lfock-1] / linear;
    emat[(nfock1-1)*20 + nfock1-1] = 0.0;
    double diag = emat[(lfock-1)*20 + lfock-1];
    if (diag < 1e-20) return;
    double c = 1.0/diag;
    for (int i=0;i<nfock;++i) for (int j=0;j<nfock;++j) emat[i*20+j]*=c;
    double evec[400];
    int l=0;
    for (int i=0;i<nfock1;++i) {
        for (int j=0;j<nfock1;++j) evec[l+j] = emat[i*20+j];
        l += nfock1;
    }
    c = diag;
    for (int i=0;i<nfock;++i) for (int j=0;j<nfock;++j) emat[i*20+j]*=c;
    double d;
    osinv(evec, nfock1, d);
    if (std::fabs(d) < 1e-6) { start = true; return; }
    if (nfock < 2) return;
    double coeffs[20];
    int il = nfock*nfock1;
    for (int i=0;i<nfock;++i) coeffs[i] = -evec[il+i];
    for (int i=0;i<linear;++i) {
        double sum = 0;
        int ii = i*mfock;
        for (int j=0;j<nfock;++j) sum += coeffs[j]*fock[j+ii];
        f[i] = sum;
    }
}
