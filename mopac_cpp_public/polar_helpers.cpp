// polar_helpers.cpp - independent helpers of polar.F90 (MOPAC2016 -> C++)
// M10 batch 11a: matrix/trace/keyword helpers + light density/work routines.
// Driver routines (polar/alphaf/betaf/beopor/non*/ng*/bmakuf/epsab) live in
// polar.cpp and require the full SCF chain (compfg->hcore->iter); that closure
// is recorded as a global gap until the SCF end-to-end batches link together.
#include "polar_helpers.h"
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fstream>

// ---- module globals (extern) ----
#include "parameters_C.h"
namespace common_arrays_C {
extern std::vector<int> nat, na, labels, nfirst, nlast;
extern std::vector<std::vector<double>> geo, coord, c;
}
namespace parameters_C {
extern double dd[N_PARAM + 1], gss[N_PARAM + 1], gsp[N_PARAM + 1], gpp[N_PARAM + 1];
extern double gp2[N_PARAM + 1], hsp[N_PARAM + 1], polvol[N_PARAM + 1];
}
namespace funcon_C { extern double ev, a0; }
namespace elemts_C { extern std::vector<std::string> elemnt; }
namespace molkst_C {
extern int numat, norbs, natoms, nclose, last;
extern std::string keywrd;
}
namespace polar_C { extern double omega; }
namespace chanel_C { extern int iw; extern std::string density_fn, pol_fn;
extern int idaf, irecln, irecst;
extern std::vector<int> ioda, ifilen; }
namespace { std::vector<std::vector<double>> w; }  // polar local work matrix
using common_arrays_C::nat; using common_arrays_C::na; using common_arrays_C::labels;
using common_arrays_C::geo; using common_arrays_C::coord; using common_arrays_C::c;
using common_arrays_C::nfirst; using common_arrays_C::nlast;
using parameters_C::polvol; using parameters_C::dd; using parameters_C::gss; using parameters_C::gsp; using parameters_C::gpp; using parameters_C::gp2; using parameters_C::hsp; using funcon_C::ev; using funcon_C::a0;
using elemts_C::elemnt; using molkst_C::numat; using molkst_C::norbs;
using molkst_C::natoms; using molkst_C::nclose; using molkst_C::last;
using molkst_C::keywrd; using polar_C::omega; using chanel_C::iw; using chanel_C::density_fn; using chanel_C::pol_fn;
using chanel_C::idaf; using chanel_C::irecln; using chanel_C::irecst;
using chanel_C::ioda; using chanel_C::ifilen;

extern double reada(const std::string&, int);
extern void mopend(const std::string&);
extern void to_screen(const std::string&);
void zerom(std::vector<std::vector<double>>& x, int m){
    for(int i=1;i<=m;++i) for(int j=1;j<=m;++j) x[i][j]=0.0;
}

void tf(std::vector<std::vector<double>>& ua, std::vector<std::vector<double>>& ga,
        std::vector<std::vector<double>>& ub, std::vector<std::vector<double>>& gb,
        std::vector<std::vector<double>>& t, int norbs){
    zerom(t, norbs);
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j){
        double sum1=0.0, sum2=0.0;
        for(int k=1;k<=norbs;++k){
            sum1 += ga[i][k]*ub[k][j] + gb[i][k]*ua[k][j] - ua[i][k]*gb[k][j] - ub[i][k]*ga[k][j];
            sum2 += ga[j][k]*ub[k][i] + gb[j][k]*ua[k][i] - ua[j][k]*gb[k][i] - ub[j][k]*ga[k][i];
        }
        t[i][j]=sum1; t[j][i]=sum2;
    }
}

void transf(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& g,
            std::vector<std::vector<double>>& c, int norb){
    std::vector<double> work(norbs+1,0.0);
    for(int i=1;i<=norb;++i){
        for(int j=1;j<=norb;++j){
            double term=0.0;
            for(int k=1;k<=norb;++k) term += f[j][k]*c[k][i];
            work[j]=term;
        }
        for(int j=1;j<=norb;++j){
            double term2=0.0;
            for(int k=1;k<=norb;++k) term2 += work[k]*c[k][j];
            g[j][i]=term2;
        }
    }
}

double trsub(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
             std::vector<std::vector<double>>& ur, int l1, int lm, int ndim){
    double sum=0.0;
    for(int i=1;i<=l1;++i) for(int k=1;k<=lm;++k){
        double suml=0.0;
        for(int l=1;l<=lm;++l) suml += x[k][l]*ur[l][i];
        sum += suml*ul[i][k];
    }
    return 2.0*sum;
}

double trudgu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim){
    double sum=0.0;
    for(int i=1;i<=l1;++i) for(int k=1;k<=lm;++k){
        double suml=0.0;
        for(int l=1;l<=lm;++l) suml += x[k][l]*ur[l][i];
        sum += suml*ul[k][i];
    }
    return 2.0*sum;
}

double trugdu(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim){
    double sum=0.0;
    for(int i=1;i<=l1;++i) for(int k=1;k<=lm;++k){
        double suml=0.0;
        for(int l=1;l<=lm;++l) suml += x[l][k]*ur[l][i];
        sum += suml*ul[i][k];
    }
    return 2.0*sum;
}

double trugud(std::vector<std::vector<double>>& ul, std::vector<std::vector<double>>& x,
              std::vector<std::vector<double>>& ur, int l1, int lm, int ndim){
    double sum=0.0;
    for(int i=1;i<=l1;++i) for(int k=1;k<=lm;++k){
        double suml=0.0;
        for(int l=1;l<=lm;++l) suml += x[k][l]*ur[i][l];
        sum += suml*ul[i][k];
    }
    return 2.0*sum;
}

double wrdkey(const std::string& keywrd_, const std::string& key, int nk,
              const std::string& refkey, int nr, double def){
    size_t i = keywrd_.find(key.substr(0,nk));
    if(i != std::string::npos){
        size_t j = keywrd_.substr(i).find(") ");
        if(j == std::string::npos) j = keywrd_.size()-i;
        size_t k = keywrd_.substr(i, j).find(refkey.substr(0,nr));
        if(k != std::string::npos) return reada(keywrd_, (int)(i+k));
        return def;
    }
    return def;
}

double aval(std::vector<std::vector<double>>& h, std::vector<std::vector<double>>& d, int norbs){
    double sum=0.0;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) sum += h[i][j]*d[j][i];
    return -sum;
}


void copym(std::vector<std::vector<double>>& h, std::vector<std::vector<double>>& f, int m){
    for(int i=1;i<=m;++i) for(int j=1;j<=m;++j) f[i][j]=h[i][j];
}

void hplusf(std::vector<std::vector<double>>& f, std::vector<std::vector<double>>& h, int norbs){
    double hartr=ev;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) f[i][j]=h[i][j]+f[i][j]/hartr;
}

void hmuf(std::vector<std::vector<double>>& h1, int id, std::vector<std::vector<double>>& coord_,
          std::vector<int>& nfirst_, std::vector<int>& nlast_, std::vector<int>& nat_, int norbs, int numat){
    zerom(h1, norbs);
    for(int i=1;i<=numat;++i){
        int ia=nfirst_[i];
        int ib=std::min(nlast_[i], ia+3);
        int ni=nat_[i];
        for(int i1=ia;i1<=ib;++i1){
            for(int j1=ia;j1<=i1;++j1){
                h1[i1][j1]=0.0;
                int io1=i1-ia, jo1=j1-ia;
                if(id==1 && jo1==0 && io1==1){ h1[i1][j1]=dd[ni]; h1[j1][i1]=dd[ni]; }
                if(id==2 && jo1==0 && io1==2){ h1[i1][j1]=dd[ni]; h1[j1][i1]=dd[ni]; }
                if(id==3 && jo1==0 && io1==3){ h1[i1][j1]=dd[ni]; h1[j1][i1]=dd[ni]; }
            }
            h1[i1][i1]=0.0;
            h1[i1][i1]+=1.8897262*coord_[id][i];
        }
    }
}
void makeuf(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& uold,
            std::vector<std::vector<double>>& g, std::vector<double>& eigs_, bool& last,
            int norbs, int nclose, double& diff, double atol){
    double hartr=ev;
    zerom(u, norbs);
    for(int k=nclose+1;k<=norbs;++k) for(int i=1;i<=nclose;++i){
        u[i][k]=hartr*g[i][k]/(eigs_[k]-eigs_[i]-omega);
        u[k][i]=hartr*g[k][i]/(eigs_[i]-eigs_[k]-omega);
    }
    diff=0.0;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) diff=std::max(std::abs(u[i][j]-uold[i][j]),diff);
    if(diff<atol) last=true;
    for(int i=1;i<=norbs;++i) for(int j=1;j<=norbs;++j) uold[i][j]=u[i][j];
}

void densf(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& c,
           std::vector<std::vector<double>>& d, std::vector<std::vector<double>>& da,
           int norbs, int nclose, std::vector<double>& w2){
    std::vector<double> w1(norbs+1,0.0);
    for(int j=1;j<=norbs;++j) for(int k=1;k<=norbs;++k){
        double sum=0.0;
        for(int l=1;l<=nclose;++l) sum+=u[k][l]*c[j][l];
        w2[(j-1)*norbs+k]=sum;
    }
    for(int i=1;i<=norbs;++i){
        for(int k=1;k<=norbs;++k){
            double sum=0.0;
            for(int l=1;l<=nclose;++l) sum+=c[i][l]*u[l][k];
            w1[k]=sum;
        }
        for(int j=1;j<=norbs;++j){
            double sum=0.0;
            for(int k=1;k<=norbs;++k) sum+=c[i][k]*w2[(j-1)*norbs+k]-w1[k]*c[j][k];
            d[i][j]=2.0*sum;
            da[i][j]=sum;
        }
    }
}

namespace parameters_C { extern double gss[N_PARAM + 1], gsp[N_PARAM + 1], gpp[N_PARAM + 1], gp2[N_PARAM + 1], hsp[N_PARAM + 1]; }
using parameters_C::gss; using parameters_C::gsp; using parameters_C::gpp;
using parameters_C::gp2; using parameters_C::hsp;

void dawrt1(std::vector<double>& v, int len, int idaf_, int ns){
    std::string dafname = "daf";
    if(!chanel_C::density_fn.empty()) dafname = chanel_C::density_fn;
    std::ofstream os(dafname, std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
    if(!os) return;
    os.seekp((long long)(ns-1)*irecln*8);
    for(int i=1;i<=len;++i) os.write((const char*)&v[i], 8);
}

void dawrit(std::vector<std::vector<double>>& v, int len, int nrec){
    int n=ioda[nrec];
    if(n<=0 || len==ifilen[nrec]){
        bool newrec=false;
        if(n<=0){
            ioda[nrec]=irecst;
            ifilen[nrec]=len;
            newrec=true;
            irecst=irecst+(len-1)/irecln+1;
            n=ioda[nrec];
        }
        int ist=(-irecln)+1;
        int ns=n;
        int lent=len;
        ist=ist+irecln;
        int ifin=ist+lent-1;
        if(ifin-ist+1>irecln) ifin=ist+irecln-1;
        int nsp=ns;
        int lenw=ifin-ist+1;
        std::vector<double> chunk(lenw+1);
        for(int k=1;k<=lenw;++k) chunk[k]=v[1][ist+k-1];
        dawrt1(chunk, lenw, idaf, nsp);
        lent-=irecln;
        ns+=1;
        n=ns;
        while(lent>=1){
            ist=ist+irecln;
            ifin=ist+lent-1;
            if(ifin-ist+1>irecln) ifin=ist+irecln-1;
            nsp=ns;
            lenw=ifin-ist+1;
            for(int k=1;k<=lenw;++k) chunk[k]=v[1][ist+k-1];
            dawrt1(chunk, lenw, idaf, nsp);
            lent-=irecln;
            ns+=1;
            n=ns;
        }
        if(newrec){
            std::ofstream os2("daf", std::ios::binary | std::ios::in | std::ios::out);
            if(os2){
                os2.seekp(0);
                os2.write((const char*)&irecst, 4);
                for(int k=1;k<=2000;++k) os2.write((const char*)&ioda[k], 4);
                for(int k=1;k<=2000;++k) os2.write((const char*)&ifilen[k], 4);
            }
        }
        return;
    }
    std::printf(" DAWRIT HAS REQUESTED A RECORD WITH LENGTH DIFFERENT THAN BEFORE - ABORT FORCED.\n DAF RECORD %5d NEW LENGTH =%5d OLD LENGTH =%5d\n", nrec, len, ifilen[nrec]);
    mopend("DAWRIT HAS REQUESTED A RECORD WITH LENGTH DIFFERENT THAN BEFORE - ABORT FORCED.");
}

void darea1(std::vector<double>& v, int len, int idaf_, int ns){
    std::string dafname = "daf";
    if(!chanel_C::density_fn.empty()) dafname = chanel_C::density_fn;
    std::ifstream is(dafname, std::ios::binary);
    if(!is) return;
    is.seekg((long long)(ns-1)*irecln*8);
    for(int i=1;i<=len;++i) is.read((char*)&v[i], 8);
}

void daread(std::vector<std::vector<double>>& v, int len, int nrec){
    int n=ioda[nrec];
    if(n!=(-1)){
        int is=(-irecln)+1;
        int ns=n;
        int lent=len;
        is=is+irecln;
        int ifin=is+lent-1;
        if(ifin-is+1>irecln) ifin=is+irecln-1;
        int nsp=ns;
        int lenw=ifin-is+1;
        std::vector<double> chunk(lenw+1);
        darea1(chunk, lenw, idaf, nsp);
        for(int k=1;k<=lenw;++k) v[1][is+k-1]=chunk[k];
        lent-=irecln;
        ns+=1;
        n=ns;
        while(lent>=1){
            is=is+irecln;
            ifin=is+lent-1;
            if(ifin-is+1>irecln) ifin=is+irecln-1;
            nsp=ns;
            lenw=ifin-is+1;
            darea1(chunk, lenw, idaf, nsp);
            for(int k=1;k<=lenw;++k) v[1][is+k-1]=chunk[k];
            lent-=irecln;
            ns+=1;
            n=ns;
        }
        return;
    }
    std::printf(" *** ERROR ***, ATTEMPT TO READ A DAF RECORD THAT WAS NEVER WRITTEN. NREC,LEN=%5d%10d\n", nrec, len);
    mopend("*** ERROR ***, ATTEMPT TO READ A DAF RECORD THAT WAS NEVER WRITTEN");
}

void openda(int irest){
    idaf=17;
    irecln=1023;
    std::string dafname = pol_fn.empty() ? std::string("daf") : pol_fn;
    std::fstream fs(dafname, std::ios::binary | std::ios::in | std::ios::out | std::ios::trunc);
    if(irest==0){
        irecst=1;
        for(int k=1;k<=2000;++k) ioda[k]=-1;
        irecst=irecst+1;
        fs.seekp(0);
        fs.write((const char*)&irecst, 4);
        for(int k=1;k<=2000;++k) fs.write((const char*)&ioda[k], 4);
        for(int k=1;k<=2000;++k) fs.write((const char*)&ifilen[k], 4);
        return;
    }
    fs.seekg(0);
    fs.read((char*)&irecst, 4);
    for(int k=1;k<=2000;++k) fs.read((char*)&ioda[k], 4);
    for(int k=1;k<=2000;++k) fs.read((char*)&ifilen[k], 4);
}

void betal1(std::vector<std::vector<double>>& u0a, std::vector<std::vector<double>>& g0a,
            std::vector<std::vector<double>>& u1b, std::vector<std::vector<double>>& g1b,
            std::vector<std::vector<double>>& u1c, std::vector<std::vector<double>>& g1c,
            int nclose, int norbs, double& term){
    double t1a=trugud(u0a,g1b,u1c,nclose,norbs,norbs);
    double t2a=trudgu(u1c,g1b,u0a,nclose,norbs,norbs);
    double t3a=trugdu(u1b,g1c,u0a,nclose,norbs,norbs);
    double t4a=trugdu(u0a,g1c,u1b,nclose,norbs,norbs);
    double t5a=trudgu(u1c,g0a,u1b,nclose,norbs,norbs);
    double t6a=trugud(u1b,g0a,u1c,nclose,norbs,norbs);
    double t1b=trugud(u0a,g1b,u1c,norbs,nclose,norbs);
    double t2b=trudgu(u1c,g1b,u0a,norbs,nclose,norbs);
    double t3b=trugdu(u1b,g1c,u0a,norbs,nclose,norbs);
    double t4b=trugdu(u0a,g1c,u1b,norbs,nclose,norbs);
    double t5b=trudgu(u1c,g0a,u1b,norbs,nclose,norbs);
    double t6b=trugud(u1b,g0a,u1c,norbs,nclose,norbs);
    term=t1b-t1a+t2b-t2a+t3a-t3b+t4a-t4b+t5b-t5a+t6b-t6a;
}

void betall(std::vector<std::vector<double>>& u2a, std::vector<std::vector<double>>& g2a,
            std::vector<std::vector<double>>& u1b, std::vector<std::vector<double>>& g1b,
            std::vector<std::vector<double>>& u1c, std::vector<std::vector<double>>& g1c,
            int nclose, int norbs, double& term){
    double t1a=trudgu(u2a,g1b,u1c,nclose,norbs,norbs);
    double t2a=trugud(u1c,g1b,u2a,nclose,norbs,norbs);
    double t3a=trugud(u1b,g1c,u2a,nclose,norbs,norbs);
    double t4a=trudgu(u2a,g1c,u1b,nclose,norbs,norbs);
    double t5a=trugdu(u1c,g2a,u1b,nclose,norbs,norbs);
    double t6a=trugdu(u1b,g2a,u1c,nclose,norbs,norbs);
    double t1b=trudgu(u2a,g1b,u1c,norbs,nclose,norbs);
    double t2b=trugud(u1c,g1b,u2a,norbs,nclose,norbs);
    double t3b=trugud(u1b,g1c,u2a,norbs,nclose,norbs);
    double t4b=trudgu(u2a,g1c,u1b,norbs,nclose,norbs);
    double t5b=trugdu(u1c,g2a,u1b,norbs,nclose,norbs);
    double t6b=trugdu(u1b,g2a,u1c,norbs,nclose,norbs);
    term=t1b-t1a+t2b-t2a+t3b-t3a+t4b-t4a+t5a-t5b+t6a-t6b;
}

void betcom(std::vector<std::vector<double>>& u1, std::vector<std::vector<double>>& g1,
            std::vector<std::vector<double>>& u2, std::vector<std::vector<double>>& g2,
            int nclose, int norbs, double& term){
    double t1a=trudgu(u2,g1,u1,nclose,norbs,norbs);
    double t2a=trugud(u1,g1,u2,nclose,norbs,norbs);
    double t3a=trugdu(u1,g2,u1,nclose,norbs,norbs);
    double t1b=trudgu(u2,g1,u1,norbs,nclose,norbs);
    double t2b=trugud(u1,g1,u2,norbs,nclose,norbs);
    double t3b=trugdu(u1,g2,u1,norbs,nclose,norbs);
    term=2.0*(t1b-t1a+t2b-t2a+t3a-t3b);
}
