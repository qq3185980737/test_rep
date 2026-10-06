// iter.cpp 鈥?SCF main loop (translated from iter.F90).
#include "iter.h"
#include "iter_C.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>

// external I-routines (declared in respective units)
extern double reada(const std::string&, int);
extern "C" { void dcopy_(int*,const double*,int*,double*,int*); }
extern void vecprt(const std::vector<double>&, int);
extern void fock2(std::vector<double>&,const std::vector<double>&,std::vector<double>&,
                  const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,
                  int,const std::vector<int>&,const std::vector<int>&,int);
extern double helect(int,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&);
extern void eigenvectors_LAPACK(std::vector<double>&,std::vector<double>&,std::vector<double>&,int);
extern double meci();
extern double second(int);
extern void timer(const char*);
extern void writmo();
extern void mopend(const char*);
extern double capcor(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&,const std::vector<double>&,const std::vector<double>&);
extern void matout(std::vector<double>&,const std::vector<double>&,int,int&,int);
extern void swap(std::vector<double>&,int,int,int,int&);
extern void density_for_GPU(std::vector<double>&,double,int,int,double,int,int,int,std::vector<double>&,int);
extern void cnvg(std::vector<double>&,std::vector<double>&,std::vector<double>&,int,double&);
extern void densit(std::vector<double>&,int,int,int,double,int,double,std::vector<double>&,int);
extern void pulay(std::vector<double>&,std::vector<double>&,int,std::vector<double>&,std::vector<double>&,std::vector<double>&,
                  int&,int&,int,bool&,double&);
extern void diag_for_GPU(std::vector<double>&,std::vector<double>&,int,std::vector<double>&,int,int);
extern void interp(int,int,int&,double,double*,double*,double*,double*,double*,double*,double*,double*,double&);
extern void memory_error(const char*);
extern void to_screen(const std::string&);
extern void chrge(const std::vector<double>&,std::vector<double>&);
extern void phase_lock(std::vector<double>&,int);
extern void delete_iter_arrays();
extern void den_in_out(int);

static void md_stop_set() {}

// globals from common_arrays_C / molkst_C (via using)
namespace common_arrays_C {
extern std::vector<double> f,fb,h,p,pa,pb,eigs,eigb,w,wk,pdiag;
extern std::vector<int> nat;
extern std::vector<int> nfirst,nlast;
extern std::vector<std::vector<double> > c;
}
// c/cb: F90 common_arrays_C holds (norbs,norbs) 2-D arrays; iter.F90 uses them as flat
// lower-half storage only via helper calls.  Use file-static flat vectors here to keep
// the SCF chain linkable; 2-D unification happens in the M03 main-chain batch.
static std::vector<double> cflat, cbflat;
namespace molkst_C {
extern int numcal,norbs,numat,nclose,nalpha,nbeta,nopen,mpack;
extern std::string keywrd;
extern bool uhf,moperr;
extern double fract,gnorm,emin,enuclr,atheat,enrgy;
extern int last,nscf,id,iscf,iflepo;
extern bool limscf;
extern int nalpha_open, nbeta_open;
}
using namespace common_arrays_C;
using namespace molkst_C;
using namespace iter_C;

static void dcopy(int n, const double* a, int, std::vector<double>& b, int bi){
    for(int i=0;i<n;++i) b[i+bi-1]=a[i];
}

// ---- DIIS (Pulay) extrapolation, enabled only via MOPAC_FAST ----
// Residual e = F*P - P*F (packed upper triangle); B_ij = <e_i,e_j>;
// solve [B 1; 1^T 0][c; lam] = [0; 1] by Gaussian elimination; F* = sum c_i F_i.
// Default (no MOPAC_FAST) path is untouched: pure official level-shift SCF.
static bool diis_solve(double* c, const std::vector<std::vector<double>>& E, int m,
                       int mpack) {
    double B[7][7];
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j) {
            double s = 0.0;
            for (int k = 1; k <= mpack; ++k) s += E[i][k] * E[j][k];
            B[i][j] = s;
        }
    double A[8][8] = {};
    double rhs[8] = {};
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) A[i][j] = B[i][j];
        A[i][m] = 1.0;
        A[m][i] = 1.0;
    }
    rhs[m] = 1.0;
    int n = m + 1;
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col + 1; r < n; ++r)
            if (std::fabs(A[r][col]) > std::fabs(A[piv][col])) piv = r;
        if (std::fabs(A[piv][col]) < 1e-12) return false;
        if (piv != col) {
            for (int j2 = col; j2 < n; ++j2) std::swap(A[col][j2], A[piv][j2]);
            std::swap(rhs[col], rhs[piv]);
        }
        for (int r = col + 1; r < n; ++r) {
            double t = A[r][col] / A[col][col];
            for (int j2 = col; j2 < n; ++j2) A[r][j2] -= t * A[col][j2];
            rhs[r] -= t * rhs[col];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        double s = rhs[r];
        for (int j2 = r + 1; j2 < n; ++j2) s -= A[r][j2] * c[j2];
        c[r] = s / A[r][r];
    }
    return true;
}

static void diis_residual(std::vector<double>& e, const std::vector<double>& f,
                          const std::vector<double>& p, int mpack, int norbs) {
#pragma omp parallel for schedule(static)
    for (int i = 1; i <= norbs; ++i) {
        for (int j = i; j <= norbs; ++j) {
            int ij = j * (j - 1) / 2 + i;
            double s = 0.0;
            for (int k = 1; k <= norbs; ++k) {
                int ik = (i <= k) ? k * (k - 1) / 2 + i : i * (i - 1) / 2 + k;
                int kj = (k <= j) ? j * (j - 1) / 2 + k : k * (k - 1) / 2 + j;
                s += f[ik] * p[kj] - p[ik] * f[kj];
            }
            e[ij] = s;
        }
    }
}

static void diis_extrapolate(std::vector<double>& f,
                             const std::vector<std::vector<double>>& Fh,
                             const double* c, int m, int mpack) {
    for (int k = 1; k <= mpack; ++k) {
        double s = 0.0;
        for (int i = 0; i < m; ++i) s += c[i] * Fh[i][k];
        f[k] = s;
    }
}

void iter(double& ee, bool fulscf, bool rand){
    auto _t0 = std::chrono::steady_clock::now();
    long long _tf2 = 0, _tdiag = 0, _tdens = 0, _tpseudo = 0, _tfull = 0;
    int _n_pseudo = 0, _n_full = 0;
    bool diis_on = getenv("MOPAC_FAST") != nullptr;
    std::vector<std::vector<double>> diis_F, diis_E;  // DIIS window (Fock, residual)
    double diis_prev_norm = 1e30;
    int diis_bad = 0;
    double selcon;
    int l=0, icalcn=0, itrmax=2000, na2el=0, na1el=0, nb1el=0, ifill=0,
         irrr=5, jalp=0, ialp=0, jbet=0, ibet=0, ihomo=1, ihomob=1, i=0, j=0,
         iemin=0, iemax=0, iredy=1, niter=0, modea=0, modeb=0;
    std::vector<double> q(numat+1,0.0), escf0(11,0.0), theta_(norbs+1,0.0);
    double plb=0, scfcrt=1e-4, pl=1.0, bshift=-80.0, pltest=0.05, trans=0.2,
           w1=0.5,w2=0.5, random=1.0, shift=0.0, shiftb=0.0, shfmax=20.0,
           ten=10.0, tenold=10.0, plchek=0.005, scorr=0.0, shfto=0.0, shftbo=0.0,
           titer0=0.0, eold=100.0, diff=0.0, enrgy_=23.060529, titer=0.0, escf=0.0,
           sellim=0.0, sum=0.0, eold_alpha_=0.0, eold_beta_=0.0;
    bool debug=false,prtfok=false,prteig=false,prtden=false,prt1el=false,minprt=false,
         newdg=false,prtpl=false,prtvec=false,camkin=false,ci=false,okpuly=false,
         oknewd=false,times=false,force=false,allcon=false,halfe=false,gs=false,
         capps=false,incitr=false,timitr=false,frst=true,bfrst=true,ready=false,
         glow=false,makea=true,makeb=true,getout=false,l_param=true;
    std::string abprt[4]={"","","ALPHA"," BETA"};
    int iopc_calcp=3;
    md_test=true; is_PARAM=false; icalcn=0; ifill=0;
    ihomo=std::max(1,nclose+nalpha); ihomob=std::max(1,nclose+nbeta);
    eold=100.0; ready=false; diff=0.0;
    if ((int)cflat.size() <= norbs * norbs) { cflat.assign(norbs * norbs + 1, 0.0); cbflat.assign(norbs * norbs + 1, 0.0); }
    if(icalcn!=numcal){
        delete_iter_arrays(); l_param=true; enrgy_=23.060529; glow=false; irrr=5; shift=0.0; icalcn=numcal; shfmax=20.0;  // F90 iter.F90:94 enrgy=fpc_9 (kcal/mol)
        debug = keywrd.find(" DEBUG")!=std::string::npos;
        minprt = (keywrd.find(" SADDLE")+0)==0;
        if(minprt) minprt = (false || debug);
        prteig = keywrd.find(" EIGS")!=std::string::npos;
        prtpl = keywrd.find(" PL ")!=std::string::npos || keywrd.find(" PLS")!=std::string::npos;
        prt1el = keywrd.find(" 1ELE")!=std::string::npos && debug;
        prtden = keywrd.find(" DENS")!=std::string::npos && debug;
        prtfok = keywrd.find(" FOCK")!=std::string::npos && debug;
        prtvec = (keywrd.find(" VEC")!=std::string::npos || keywrd.find(" ALLVEC")!=std::string::npos) && debug;
        debug = keywrd.find(" ITER")!=std::string::npos;
        newdg=false; camkin=false; plchek=0.005; pl=1.0; plb=0; bshift=-80.0; shift=1.0;  // F90 iter.F90:97 shift=0 then :121 shift=1 -> final 1.0
        if (diis_on) shift = 0.0;  // DIIS needs the pure Fock (no level shift) for extrapolation
        shfto=0; shftbo=0; itrmax=2000; na2el=nclose; na1el=nalpha+nopen; nb1el=nbeta+nopen;
        if(keywrd.find(" FILL")!=std::string::npos) ifill=-(int)reada(keywrd,(int)keywrd.find(" FILL"));
        if(keywrd.find(" SHIFT")!=std::string::npos) bshift=-reada(keywrd,(int)keywrd.find(" SHIFT"));
        if(std::abs(bshift)>1e-20) ten=bshift;
        if(keywrd.find(" ITRY")!=std::string::npos) itrmax=(int)reada(keywrd,(int)keywrd.find(" ITRY"));
        ci = keywrd.find(" MICROS")!=std::string::npos || keywrd.find(" C.I.")!=std::string::npos;
        okpuly = keywrd.find(" PULAY")!=std::string::npos;
        oknewd = std::abs(bshift)<0.001;
        if(camkin && std::abs(bshift)>1e-5) bshift=4.44;
        times = keywrd.find(" TIMES")!=std::string::npos; timitr=times;
        force = keywrd.find(" FORCE")!=std::string::npos;
        gs = (keywrd.find(" TS")==std::string::npos && keywrd.find(" NLLSQ")==std::string::npos && keywrd.find(" SIGMA")==std::string::npos) && !force;
        allcon = okpuly || camkin;
        j=0;
        for(int ii=1;ii<=numat;++ii) if(nat[ii]==102) ++j;
        capps=j>0; iscf=1; trans=0.2;
        if(keywrd.find(" OLDENS")!=std::string::npos){
            den_in_out(0); if(moperr) return;
            if(uhf){ for(int ii=1;ii<=mpack;++ii) pold[ii]=pa[ii]; for(int ii=1;ii<=mpack;++ii) pbold[ii]=pb[ii]; }
            else { for(int ii=1;ii<=mpack;++ii) pold[ii]=pa[ii]*2.0; }
        } else {
            if(!is_PARAM){
                for(int ii=1;ii<=mpack;++ii){ p[ii]=0; pa[ii]=0; pb[ii]=0; }
                w1=na1el/(na1el+1e-6+nb1el); w2=1.0-w1;
                if(w1<1e-6) w1=0.5; if(w2<1e-6) w2=0.5;
                random=1.0;
                glow = glow || (gnorm<2.0 && gnorm>1e-9);
                if(!glow && uhf && na1el==nb1el) random=1.1;
                for(int ii=1;ii<=norbs;++ii){
                    int jj=(ii*(ii+1))/2; p[jj]=pdiag[ii]; // F90 iter.F90:181 p(j)=pdiag(i)
                    pa[jj]=p[jj]*w1*random; random=1.0/random; pb[jj]=p[jj]*w2*random;
                }
                if(uhf){ for(int ii=1;ii<=norbs;++ii){ random=1.0/random; int jj=(ii*(ii+1))/2; pb[jj]=p[jj]*w2*random; } }
            }
            for(int ii=1;ii<=mpack;++ii) pold[ii]=pa[ii];
            if(uhf){ for(int ii=1;ii<=mpack;++ii) pbold[ii]=pb[ii]; }
            for(int ii=1;ii<=norbs;++ii) pold2[ii]=pold[(ii*(ii+1))/2];
        }
        halfe = (nopen!=nclose && std::abs(fract-2.0)>1e-20 && std::abs(fract)>1e-20);
        if (diis_on && (uhf || halfe || ci)) diis_on = false;  // DIIS: RHF non-CI only
        if(halfe){ iopc_calcp=3; if(lgpu) iopc_calcp=2; }
        else { iopc_calcp=5; if(lgpu) iopc_calcp=4; }
        if(gs) gs = !halfe && !ci;
        scfcrt=1e-4;
        if((keywrd.find(" NLLSQ")!=std::string::npos || keywrd.find(" SIGMA")!=std::string::npos || keywrd.find(" TS")!=std::string::npos) || force) scfcrt*=0.001;
        else if(keywrd.find(" PRECISE")!=std::string::npos || nopen!=nclose) scfcrt*=0.01;
        if(keywrd.find(" POLAR")!=std::string::npos) scfcrt=std::min(1e-6,scfcrt);
        scfcrt=std::max(scfcrt,1e-12);
        int ii=keywrd.find(" SCFCRT"), jj=keywrd.find(" RELSCF");
        if(ii!=-1) scfcrt=reada(keywrd,ii);
        else if(jj!=-1) scfcrt=reada(keywrd,jj)*scfcrt;
        if(id==3) scfcrt=scfcrt*numat/20;
    } else if(nscf>0 && !uhf){
        for(int ii=1;ii<=mpack;++ii) pb[ii]=pa.data()[ii];;
        for(int ii=1;ii<=mpack;++ii) p[ii]=2.0*pa[ii];
    }
    makea=true; makeb=true; iemin=0; iemax=0;
    if(irrr!=5){
        if(uhf){
            for(int ii=1;ii<=mpack;++ii) pold[ii]=pa.data()[ii];; for(int ii=1;ii<=mpack;++ii) pbold[ii]=pb.data()[ii];;
            for(int ii=1;ii<=norbs;++ii){ pold2[ii]=pa[(ii*(ii+1))/2]; pbold2[ii]=pb[(ii*(ii+1))/2]; }
        } else {
            for(int ii=1;ii<=mpack;++ii) pold[ii]=p.data()[ii];;
            for(int ii=1;ii<=norbs;++ii) pold2[ii]=p[(ii*(ii+1))/2];
        }
    }
    camkin = keywrd.find(" KING")!=std::string::npos || keywrd.find(" CAMP")!=std::string::npos;
    if(!fulscf) shift=0;
    if(newdg) newdg=std::abs(bshift)<0.001;
    if(last==1) newdg=false;
    selcon=scfcrt;
    if(gs) selcon=std::min(std::min(scfcrt*100.0,0.1),std::max(scfcrt*std::pow(gnorm,3)*std::pow(10.0,-id*3),scfcrt));
    pltest=0.05*std::sqrt(std::abs(selcon));
    if(nalpha!=nbeta && uhf) pltest=0.001;
    if(prt1el) vecprt(h,norbs);
    iredy=1;
L180:
    niter=0; frst=true;
    if (diis_on) { diis_F.clear(); diis_E.clear(); diis_prev_norm = 1e30; diis_bad = 0; }
    modea = camkin?1:0; modeb = camkin?1:0;
    bfrst=true;
L250:
    incitr = (modea!=3 && modeb!=3);
    if(incitr) ++niter;
    if(timitr){ titer=second(1); titer0=titer; }
    if(std::abs(shift)>1e-10 && bshift!=0.0){
        l=0;
        if(niter>1){
            if(newdg && !(halfe||camkin)){
                if(diff>0){ shift=1.0; if(diff>1) newdg=false; } else shift=-0.1;
            } else shift=ten+eigs[ihomo+1]-eigs[ihomo]+shift;
            if(diff>0){
                if(shift>4.0) shfmax=4.5;
                if(shift>shfmax) shfmax=std::max(shfmax-0.5,0.0);
            }
            shift=std::max(-20.0,std::min(shfmax,shift));
            if(std::abs(shift-shfmax)<1e-5) shfmax+=0.01;
            if(okpuly || std::abs(bshift-4.44)<1e-5){ shift=-8.0; if(newdg) shift=0; }
            if(uhf){
                if(newdg && !(halfe||camkin)) shiftb=ten-tenold;
                else shiftb=ten+eigs[ihomob+1]-eigs[ihomob]+shiftb;
                if(diff>0) shiftb=std::min(4.0,shiftb);
                shiftb=std::max(-20.0,std::min(shfmax,shiftb));
                if(okpuly || std::abs(bshift-4.44)<1e-5){ shiftb=-8.0; if(newdg) shiftb=0; }
                for(int ii=ihomob+1;ii<=norbs;++ii) eigb[ii]+=shiftb;
            }
        }
        tenold=ten;
        if(pl>plchek){ shftbo=shiftb; shfto=shift; } else { shiftb=shftbo; shift=shfto; }
        if(id!=0) shift=0;
        for(int ii=ihomo+1;ii<=norbs;++ii) eigs[ii]+=shift;
        if(id!=0) shift=-80;
        for(int ii=1;ii<=mpack;++ii) f[ii]=h[ii]+shift*pa[ii]+1e-16*ii;
        for(int ii=1;ii<=norbs;++ii) f[(ii*(ii+1))/2]-=shift;
    } else if(last==0 && niter<2 && fulscf){
        random=0.001; glow = glow || (gnorm<2.0 && gnorm>1e-9);
        if(glow) random=0;
        for(int ii=1;ii<=mpack;++ii){ random=-random; f[ii]=h[ii]+random; }
    } else { for(int ii=1;ii<=mpack;++ii) f[ii]=h.data()[ii];; }
L320:
    {
        auto _b0 = std::chrono::steady_clock::now();
        if(id!=0) fock2(f,p,pa,w,w,wk,numat,nfirst,nlast,2);
        else      fock2(f,p,pa,w,w,w,numat,nfirst,nlast,2);
        _tf2 += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - _b0).count();
    }
    if (diis_on && !getout) {
        // DIIS residual uses the OLD density (the one F_n was built from):
        //   e_n = F_n*P_{n-1} - P_{n-1}*F_n   (P_n would make e ~ 0 by construction)
        std::vector<double> e(mpack + 1, 0.0);
        diis_residual(e, f, p, mpack, norbs);
        double norm = 0.0;
        for (int k = 1; k <= mpack; ++k) norm += e[k] * e[k];
        norm = std::sqrt(norm);
        fprintf(stderr, "[DIISe] it=%d |e|=%.6e prev=%.6e bad=%d\n", niter, norm,
                diis_prev_norm, diis_bad);
        if (norm > diis_prev_norm * 1.2) ++diis_bad; else diis_bad = 0;
        if (diis_bad >= 3) { diis_F.clear(); diis_E.clear(); diis_bad = 0; }
        diis_prev_norm = norm;
        if (diis_F.size() >= 6) { diis_F.erase(diis_F.begin()); diis_E.erase(diis_E.begin()); }
        diis_F.push_back(f);
        diis_E.push_back(e);
        // residual-based convergence: F*P-P*F small => self-consistent
        if (norm < 1e-3 && niter > 5) { incitr = true; goto L470; }
        // extrapolate F* = sum c_i F_i, damped: f = 0.7*F* + 0.3*F_n (stable).
        // Delayed start (niter>=5) and |c_i|<=2 guard against wild early extrapolation.
        if (diis_F.size() >= 2 && niter >= 5) {
            double c[8];
            int m = (int)diis_F.size();
            bool ok = diis_solve(c, diis_E, m, mpack);
            if (ok) for (int i = 0; i < m; ++i) if (std::fabs(c[i]) > 2.0) ok = false;
            if (!ok && m > 2) {
                diis_F.erase(diis_F.begin());
                diis_E.erase(diis_E.begin());
                m = (int)diis_F.size();
                ok = (m >= 2) && diis_solve(c, diis_E, m, mpack);
                if (ok) for (int i = 0; i < m; ++i) if (std::fabs(c[i]) > 2.0) ok = false;
            }
            if (ok && m >= 2) {
                fprintf(stderr, "[DIIS] it=%d m=%d c=", niter, m);
                for (int i = 0; i < m; ++i) fprintf(stderr, "%.4f ", c[i]);
                fprintf(stderr, "\n");
                std::vector<double> fstar(mpack + 1);
                diis_extrapolate(fstar, diis_F, c, m, mpack);
                for (int k = 1; k <= mpack; ++k) f[k] = 0.7 * fstar[k] + 0.3 * f[k];
            }
        }
    }
    if(prtfok) vecprt(f,norbs);
    if(uhf){
        if(shiftb!=0.0){
            l=0;
            for(int ii=1;ii<=norbs;++ii){
                for(int k=1;k<=ii;++k) fb[l+k]=h[l+k]+shiftb*pb[l+k];
                l+=ii; fb[l]-=shiftb;
            }
        } else if(rand && last==0 && niter<2 && fulscf){
            random=0.001; if(glow) random=0;
            for(int ii=1;ii<=mpack;++ii){ random=-random; fb[ii]=h[ii]+random; }
        } else for(int ii=1;ii<=mpack;++ii) fb[ii]=h.data()[ii];;
        auto _b1 = std::chrono::steady_clock::now();
        if(id!=0) fock2(fb,p,pb,w,w,wk,numat,nfirst,nlast,2);
        else      fock2(fb,p,pb,w,w,w,numat,nfirst,nlast,2);
        _tf2 += std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::steady_clock::now() - _b1).count();
        if(prtfok) vecprt(fb,norbs);
    }
    if(!fulscf) goto L600;
    if(irrr==0){ for(int ii=1;ii<=norbs;++ii) f[(ii*(ii+1))/2]*=0.5; }
    irrr=2;
    if(niter>=itrmax){
        if(diff<1e-3 && pl<1e-4 && !force){ incitr=true; getout=true; goto L410; }
        if(md_test){
            okpuly=true; camkin=!halfe; allcon=true; bshift=4.44;
            iredy=-4; eold=100.0; newdg=false; md_test=false; goto L180;
        }
        iscf=2; md_stop_set();
        writmo(); mopend("UNABLE TO ACHIEVE SELF-CONSISTENCE"); return;
    }
    ee=helect(norbs,pa,h,f);
    if(uhf) ee+=helect(norbs,pb,h,fb); else ee*=2.0;
    if(capps) ee+=capcor(nat,nfirst,nlast,p,h);
    if(uhf){
        if(bshift!=0.0){
            if(nalpha_open>nalpha) scorr=shift*(nalpha_open-nalpha)*enrgy_*0.5*fract*(1.0-fract);
            else scorr=shift*(nbeta_open-nbeta)*enrgy_*0.5*fract*(1.0-fract);
        }
    } else if(bshift!=0.0) scorr=shift*(nopen-nclose)*enrgy_*0.25*fract*(2.0-fract);
    escf=(ee+enuclr)*enrgy_+atheat+scorr;
    getout=false;
L410:
    if(incitr){
        if(getout) goto L470;
        diff=escf-eold;
        if(diff>0) ten-=1.0; else ten=ten*0.975+0.05;
        sellim=std::max(selcon,1e-15*std::max(std::abs(ee),1.0));
        bool conv_nw = (niter>4 && (pl==0.0 || (pl<pltest && std::abs(diff)<sellim)) && ready);
        if(conv_nw){ goto L470; }
        goto L490;
L470:
        if(std::abs(shift)<1e-10) goto L600;
        shift=0; shiftb=0;
        for(int ii=1;ii<=mpack;++ii) f[ii]=h[ii];
        makea=true; makeb=true; goto L320;
L490:
        if(limscf && emin!=0.0 && !(ci||halfe)){
            if(escf<emin){
                iemax=0; iemin=std::min(5,iemin+1);
                for(int ii=1;ii<iemin;++ii) escf0[ii]=escf0[ii+1];
                escf0[iemin]=escf;
                if(iemin>3){
                    for(int ii=2;ii<=iemin;++ii)
                        if(std::abs(escf0[ii]-escf0[ii-1])>0.05*(emin-escf)) goto L540;
                    incitr=true; getout=true; goto L410;
                }
            } else {
                iemin=0; iemax=std::min(5,iemax+1);
                for(int ii=1;ii<iemax;++ii) escf0[ii]=escf0[ii+1];
                escf0[iemax]=escf;
                if(iemax>3){
                    for(int ii=2;ii<=iemax;++ii)
                        if(std::abs(escf0[ii]-escf0[ii-1])>0.05*(escf-emin)) goto L540;
                    incitr=true; getout=true; goto L410;
                }
            }
        }
L540:
        ready = iredy>0 && (std::abs(diff)<sellim*10.0 || pl==0.0);
        ++iredy;
    }
    if(incitr) eold=escf;
    if(niter>2 && camkin && makea){
        if(vec_ai.empty()){
            vec_ai.assign((norbs+1)*(norbs+1),0.0); fock_ai.assign((norbs+1)*(norbs+1),0.0);
            p_ai.assign((norbs+1)*(norbs+1),0.0); h_ai.assign(norbs*norbs+1,0.0); vecl_ai.assign(norbs*norbs+1,0.0);
        }
        interp(na1el,norbs-na1el,modea,escf/enrgy_,f.data(),cflat.data(),theta_.data(),vec_ai.data(),fock_ai.data(),p_ai.data(),h_ai.data(),vecl_ai.data(),eold_alpha_);
    }
    makeb=false;
    if(modea!=3){
        makeb=true;
        if(nscf==2){}
        if(newdg){
            if(okpuly && makea && iredy>1)
                pulay(f,pa,norbs,pold,pold2,pold3,jalp,ialp,6*mpack,frst,pl);
            if(halfe||camkin) eigenvectors_LAPACK(cflat,f,eigs,norbs);
            else { auto _b2 = std::chrono::steady_clock::now();
            diag_for_GPU(f,cflat,na1el,eigs,norbs,mpack);
            _tdiag += std::chrono::duration_cast<std::chrono::microseconds>(
                          std::chrono::steady_clock::now() - _b2).count();
            ++_n_pseudo;
            _tpseudo += std::chrono::duration_cast<std::chrono::microseconds>(
                          std::chrono::steady_clock::now() - _b2).count();
            }
        } else {
            auto _b3 = std::chrono::steady_clock::now();
            eigenvectors_LAPACK(cflat,f,eigs,norbs);
            if (niter == 1 && std::getenv("EIGDBG"))
                std::fprintf(stderr, "[EIGDBG1] %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
                             eigs[1], eigs[2], eigs[3], eigs[4], eigs[5], eigs[6], eigs[7], eigs[8]);
            if ((niter == 2 || niter == 3 || niter == 4 || niter == 5 || niter == 10 || niter == 20 || niter == 30) && std::getenv("EIGDBG"))
                std::fprintf(stderr, "[EIGDBG%d] %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n",
                             niter, eigs[1], eigs[2], eigs[3], eigs[4], eigs[5], eigs[6], eigs[7], eigs[8]);
            _tdiag += std::chrono::duration_cast<std::chrono::microseconds>(
                          std::chrono::steady_clock::now() - _b3).count();
            ++_n_full;
            _tfull += std::chrono::duration_cast<std::chrono::microseconds>(
                          std::chrono::steady_clock::now() - _b3).count();
        }
        if (std::getenv("SCFPROBE")) {
            std::printf(" SCFPROBE iter%6d %d\n", niter, (int)newdg);
            std::printf("  ei ");
            for (int ii = 1; ii <= 8 && ii <= norbs; ++ii) std::printf("%13.7f", eigs[ii]);
            std::printf("\n");
        }
        if(prtvec){ j=uhf?2:1; matout(cflat,eigs,norbs,norbs,norbs); }
    }
    if(ifill!=0) swap(cflat,norbs,norbs,na2el,ifill);
    if(uhf){
        density_for_GPU(cflat,fract,nalpha,nalpha_open,1.0,mpack,norbs,1,pa,iopc_calcp);
        if(modea!=3 && !(newdg && okpuly)){ i=niter; if(camkin) i=7; cnvg(pa,pold,pold2,i,pl); }
    } else {
        auto _b4 = std::chrono::steady_clock::now();
        if(halfe) densit(cflat,norbs,norbs,na2el,2.0,na1el,fract,p,1);
        else density_for_GPU(cflat,fract,na2el,na1el,2.0,mpack,norbs,1,p,iopc_calcp);
        _tdens += std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::steady_clock::now() - _b4).count();
        if(modea!=3 && !(newdg && okpuly)) cnvg(p,pold,pold2,niter,pl);
    }
    if(uhf){
        if(niter>2 && camkin && makeb){
            if(vec_bi.empty()){
                vec_bi.assign((norbs+1)*(norbs+1),0.0); fock_bi.assign((norbs+1)*(norbs+1),0.0);
                p_bi.assign((norbs+1)*(norbs+1),0.0); h_bi.assign(norbs*norbs+1,0.0); vecl_bi.assign(norbs*norbs+1,0.0);
            }
            interp(nb1el,norbs-nb1el,modeb,escf/enrgy_,fb.data(),cbflat.data(),theta_.data(),vec_bi.data(),fock_bi.data(),p_bi.data(),h_bi.data(),vecl_bi.data(),eold_beta_);
        }
        makea=false;
        if(modeb!=3){
            makea=true;
            if(newdg){
                if(okpuly && makeb && iredy>1)
                    pulay(fb,pb,norbs,pbold,pbold2,pbold3,jbet,ibet,6*mpack,bfrst,plb);
                if(halfe||camkin) eigenvectors_LAPACK(cbflat,fb,eigb,norbs);
                else diag_for_GPU(fb,cbflat,nb1el,eigb,norbs,mpack);
            } else eigenvectors_LAPACK(cbflat,fb,eigb,norbs);
            if(prtvec) matout(cbflat,eigb,norbs,norbs,norbs);
        }
        density_for_GPU(cbflat,fract,nbeta,nbeta_open,1.0,mpack,norbs,1,pb,iopc_calcp);
        if(!(newdg && okpuly)){ i=niter; if(camkin) i=7; cnvg(pb,pbold,pbold2,i,plb); }
    }
    if(uhf){ for(int ii=1;ii<=mpack;++ii) p[ii]=pa[ii]+pb[ii]; }
    else { for(int ii=1;ii<=mpack;++ii){ pa[ii]=p[ii]*0.5; pb[ii]=pa[ii]; } }
    if(itrmax<3) return;
    if (std::getenv("EIGDBG"))
        std::fprintf(stderr, "[SCFDBG] niter=%d pl=%.3e diff=%.3e ready=%d newdg=%d shift=%.3e sellim=%.3e escf=%.9f\n",
                     niter, pl, diff, (int)ready, (int)newdg, shift, sellim, escf);
    oknewd = pl<sellim || oknewd;
    newdg = (pl<trans && oknewd) || newdg;
    if(pl<trans*0.3333) oknewd=true;
    goto L250;
L600:
    ee=helect(norbs,pa,h,f);
    if(uhf) ee+=helect(norbs,pb,h,fb); else ee*=2.0;
    if(capps) ee+=capcor(nat,nfirst,nlast,p,h);
    if(nscf==0 || last==1 || ci || halfe){
        for(int ii=1;ii<=mpack;++ii) pold[ii]=f.data()[ii];;
        auto _b5 = std::chrono::steady_clock::now();
        eigenvectors_LAPACK(cflat,pold,eigs,norbs);
        _tdiag += std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::steady_clock::now() - _b5).count();
        // Sync flat column-major c (0-based physical) into common 2-D c
        // for meci/halfe, which read c[row][col] per official Fortran.
        if ((int)common_arrays_C::c.size() < norbs + 1)
            common_arrays_C::c.assign(norbs + 1,
                                      std::vector<double>(norbs + 1, 0.0));
        for (int col = 1; col <= norbs; ++col)
            for (int row = 1; row <= norbs; ++row)
                common_arrays_C::c[row][col] = cflat[(col - 1) * norbs + (row - 1)];
        if(last==1) phase_lock(cflat,norbs);
        if(uhf){
            for(int ii=1;ii<=mpack;++ii) pold[ii]=fb.data()[ii];; eigenvectors_LAPACK(cbflat,pold,eigb,norbs);
            if(last==1) phase_lock(cbflat,norbs);
            for(int ii=1;ii<=mpack;++ii) pold[ii]=pa.data()[ii];;
        } else for(int ii=1;ii<=mpack;++ii) pold[ii]=p.data()[ii];;
        if(ci||halfe){
            sum=meci();
            if(moperr) return;
            ee+=sum;
        }
    }
    ++nscf;
    if(allcon && std::abs(bshift-4.44)<1e-7){
        camkin=false; allcon=false; newdg=false; bshift=-10; okpuly=false;
    }
    shift=1.0;
    escf=(ee+enuclr)*enrgy_+atheat;
    if(emin==0.0) emin=escf; else emin=std::min(emin,escf);
    fprintf(stderr, "[ITPROF] fock2=%.1fms diag=%.1fms(pseudo:%d=%.1fms full:%d=%.1fms) dens=%.1fms total=%.1fms\n",
            _tf2 / 1e3, _tdiag / 1e3, _n_pseudo, _tpseudo / 1e3, _n_full, _tfull / 1e3,
            _tdens / 1e3,
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now() - _t0).count() / 1e3);
}

// delete_iter_arrays: F90 deallocates scratch arrays; C++ RAII makes this a no-op.
void delete_iter_arrays() {}
