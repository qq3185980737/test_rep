// selmos.cpp
#include "selmos.h"
#include "MOZYME_C.h"
#include <cstdio>
#include <algorithm>
using namespace MOZYME_C;
extern int iw;

void selmos(int nmos_loc, int* nc, int* ic, int n01, double* c, int n02,
            int* nnc_loc, int* ncmo, int ln, int mn, int* iws, int* jopt,
            int* ncnew, int* ncmnew, int* nncnew, int mode) {
    int l = n02 - mn;
    for (int i=mn;i>=1;--i) c[i+l]=c[i];
    for (int i=1;i<=nmos_loc;++i) ncmo[i]+=l;
    l = n01 - ln;
    for (int i=ln;i>=1;--i) ic[i+l]=ic[i];
    for (int i=1;i<=nmos_loc;++i) nnc_loc[i]+=l;
    int jbot=0, nbot=0, jtop=n02-mn-1, ntop=n01-ln-1;
    int ibot=0, itop=nmos_loc+1, nreal=0;
    for (int i=nmos_loc;i>=2;--i) iws[i]-=iws[i-1];
    for (int i=1;i<=nmos_loc;++i) {
        int ll=nnc_loc[i], l1=ic[1+ll], l2=ic[2+ll];
        bool sel=false;
        for (int m=1;m<=numred;++m)
            if (l1==jopt[m]||l2==jopt[m]) { sel=true; break; }
        if (!sel) {
            --itop;
            nncnew[itop]=nnc_loc[i]; ncnew[itop]=nc[i]; ncmnew[itop]=ncmo[i];
            continue;
        }
        if (nbot+nc[i]>ntop || jbot+iws[i]>jtop) {
            // compct stub
        }
        ++ibot;
        nncnew[ibot]=nbot; ncnew[ibot]=nc[i]; ncmnew[ibot]=jbot;
        for (int n=1;n<=nc[i];++n) ic[nbot+n]=ic[ll+n];
        nbot+=nc[i]; nc[i]=0;
        int j=ncmo[i]; int ncoefs=iws[i];
        for (int n=1;n<=ncoefs;++n) c[jbot+n]=c[j+n];
        jbot+=ncoefs; ++nreal;
    }
    for (int i=1;i<=ibot;++i){nc[i]=ncnew[i];nnc_loc[i]=nncnew[i];ncmo[i]=ncmnew[i];}
    int j=nmos_loc+1;
    for (int i=itop;i<=nmos_loc;++i){ --j; nc[i]=ncnew[j];nnc_loc[i]=nncnew[j];ncmo[i]=ncmnew[j];}
    if (jbot>jtop) { std::printf("BUG IN SELMOS\n"); }
    if (nbot>ntop) { std::printf("BUG IN SELMOS\n"); }
    if (mode==1) nelred=2*nreal; else norred=nreal+nelred/2;
}
