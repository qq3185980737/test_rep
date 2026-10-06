// tidy.cpp
#include "tidy.h"
#include "MOZYME_C.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
using namespace MOZYME_C;
using namespace molkst_C;
using namespace chanel_C;
extern int iw;
extern void selmos(int, int*, int*, int, double*, int, int*, int*, int, int,
                   int*, int*, int*, int*, int*, int);
extern void memory_error(const char*);
extern void pinout(int);

void tidy(int nmos_loc, int* nc, int* ic, int& n01, double* c, int n02,
          int* nnc_loc, int* ncmo, int& ln, int& mn, int mode) {
    static int icalcn=0, ireset=0;
    static bool debug=false, large=false;
    static int imode[3]={0,0,0};
    std::vector<int> iused(norbs+1,0), jused(norbs+1,0), kused(norbs+1,0), lused(norbs+1,0);
    if (icalcn!=numcal) {
        icalcn=numcal;
        debug = keywrd.find(" TIDY") != std::string::npos;
        large = keywrd.find(" LARGE") != std::string::npos;
    }
    if (nmos_loc==0) return;
    int isnew=0, nmol=0;
    if (nmol!=numcal) { imode[1]=step_num; imode[2]=step_num; ireset=0; }
    int momax=0, momin=10000;
    while (true) {
        momax=0; --ireset;
        if (ireset==0) thresh*=100.0;
        momin=10000; ln=0; mn=0; int in=0, li=0;
        for (int i=1;i<=nmos_loc;++i) {
            int ll=nnc_loc[i]; li-=ll; int mm=ncmo[i];
            ncmo[i]=mn; nnc_loc[i]=0;
            if (i!=1) nnc_loc[i]=nnc_loc[i-1]+in;
            in=0;
            momax=std::max(momax,nc[i]);
            momin=std::min(momin,nc[i]);
            for (int j1=1;j1<=nc[i];++j1) {
                ++ll; double sum=0; int j=ic[ll];
                for (int k=1;k<=iorbs[j];++k){++mm; sum+=c[mm]*c[mm];}
                if (sum>thresh) {
                    mm-=iorbs[j];
                    for (int k=1;k<=iorbs[j];++k){++mn;++mm;c[mn]=c[mm];}
                    ++in; ++ln; ic[ln]=j;
                }
            }
            nc[i]=in; iused[i]=mn; li+=ll;
        }
        if (isnew!=0||step_num<=1||step_num==imode[mode]) break;
        isnew=2;
        if (numred>=numat-1) isnew=1;
        imode[mode]=step_num;
        if (isnew==2) selmos(nmos_loc,nc,ic,n01,c,n02,nnc_loc,ncmo,ln,mn,iused.data(),jopt.data(),jused.data(),kused.data(),lused.data(),mode);
        else break;
    }
    int ispace=(n01-ln)/nmos_loc;
    int jspace=(n02-mn)/nmos_loc;
    int i=n01/nmos_loc;
    if (i+ispace<numat && ispace<std::max(i/5,20)) { pinout(1); moperr=true; return; }
    if (ispace>numat && jspace>nmos_loc) { ispace=0; jspace=0; }
    int jtop=n01, mtop=n02;
    int jsav=n01-ln-ispace*nmos_loc, msav=n02-mn-jspace*nmos_loc;
    for (int ii=nmos_loc;ii>=2;--ii) {
        int j=std::min(nc[ii]+ispace,numat);
        if (j<nc[ii]+ispace) { jsav+=nc[ii]+ispace-j; jtop-=j; }
        else { int jdash=std::min(numat,std::max(jsav,0)+j); jtop-=jdash; jsav=jsav-jdash+j; }
        int l=nnc_loc[ii]; nnc_loc[ii]=jtop;
        for (int k=nc[ii];k>=1;--k) ic[jtop+k]=ic[l+k];
        int n=iused[ii]-iused[ii-1];
        j=std::min(n+jspace,norbs);
        if (j<n+jspace) { msav+=n+jspace-j; mtop-=j; }
        else { int mdash=std::min(norbs,std::max(msav,0)+j); mtop-=mdash; msav=msav-mdash+j; }
        ncmo[ii]=mtop;
        for (int k=n;k>=1;--k) c[k+mtop]=c[k+iused[ii-1]];
    }
    if (mode==2) nmol=numcal;
}
