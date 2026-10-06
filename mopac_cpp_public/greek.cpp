// greek.cpp
#include "greek.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
using common_arrays_C::txtatm;
using common_arrays_C::nat;
using common_arrays_C::nbonds;
using common_arrays_C::ibonds;
void txtype(int nnext, int* inext, char type);

void greek(int n1) {
    int ihcr = MOZYME_C::nbackb[0];
    int ico  = MOZYME_C::nbackb[1];
    int jatom1 = 0;
    (void)MOZYME_C::nbackb[3];
    int i,j;
    for (i=1;i<=nbonds[ihcr];++i)
        if (nat[ibonds[i][ihcr]]==8) goto done;
    txtatm[ihcr][14]='A';
done:
    int beta=0;
    for (i=1;i<=nbonds[ihcr];++i) {
        j=ibonds[i][ihcr];
        if (nat[j]==6 && j!=ico) { beta=j; break; }
    }
    if (!beta) return;
    txtatm[beta][14]='B';
    bool Trp=txtatm[beta].substr(17,3)=="TRP";
    bool Ile=txtatm[beta].substr(17,3)=="ILE";
    std::vector<int> icurr(81),inext(81),iprev(81);
    icurr[1]=beta; int ncurr=1; iprev[1]=ihcr;
    char types[23]={' ','G','D','E','Z','H','T','I','K','L','M',
                    'N','X','O','P','R','S','T','U','F','C','Y','W'};
    for (int jjj=1;jjj<=22;++jjj) {
        int nprev=ncurr; int nnext=0;
        for (int ii=1;ii<=ncurr;++ii) {
            j=icurr[ii];
            if (j==n1) continue;
            for (int l=1;l<=nbonds[j];++l) {
                int m=ibonds[l][j];
                bool skip=false;
                for (int ll=1;ll<=nprev;++ll) if (m==iprev[ll]) { skip=true; break; }
                if (skip) continue;
                for (int ll=1;ll<=ncurr;++ll) if (m==icurr[ll]) { skip=true; break; }
                if (skip) continue;
                if (m!=jatom1) { ++nnext; inext[nnext]=m; }
                if (nat[m]==7 && txtatm[ihcr].substr(17,3)=="PRO")
                    nnext = std::max(0,nnext-1);
            }
        }
        if (nnext==0) break;
        if (Trp && jjj==2) {
            j=inext[1];
            int l; for (l=1;l<=nbonds[j];++l)
                if (nat[ibonds[l][j]]==7) break;
            if (l>nbonds[j]) {
                int tmp=inext[1]; inext[1]=inext[2]; inext[2]=tmp;
            }
        } else if (Ile && jjj==1) {
            j=inext[1]; int m=0;
            for (int l=1;l<=nbonds[j];++l)
                if (nat[ibonds[l][j]]==6) ++m;
            if (m==1) { int tmp=inext[1]; inext[1]=inext[2]; inext[2]=tmp; }
        }
        txtype(nnext, inext.data(), types[jjj]);
        for (int k=1;k<=ncurr;++k) iprev[k]=icurr[k];
        for (int k=1;k<=nnext;++k) icurr[k]=inext[k];
        ncurr=nnext;
    }
}
void txtype(int, int*, char) {}
