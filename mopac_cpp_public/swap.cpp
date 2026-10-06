// swap.cpp
#include "swap.h"
#include <cmath>
#include <algorithm>
namespace iter_C { std::vector<double> psi, stdpsi; }
using iter_C::psi;
using iter_C::stdpsi;
static double ddot(int n,const double*a,const double*b){double s=0;for(int i=0;i<n;++i)s+=a[i]*b[i];return s;}
extern int iw;
#include <cstdio>

void swap(std::vector<double>& c, int n, int mdim, int nocc, int& ifill) {
    swap(c.data(), n, mdim, nocc, ifill);
}
void swap(double* c, int n, int mdim, int nocc, int& ifill) {
    if (ifill <= 0) {
        psi.assign(n, 0.0);
        stdpsi.assign(n, 0.0);
        ifill = -ifill;
        for (int i=1;i<=n;++i){ stdpsi[i-1]=c[(ifill-1)*mdim+(i-1)]; psi[i-1]=c[(ifill-1)*mdim+(i-1)]; }
        return;
    }
    double sum = ddot(n, psi.data(), &c[(ifill-1)*mdim]);
    if (std::fabs(sum) <= 0.707106781187) {
        double summax=0; int jfill=1;
        for (int f=1;f<=n;++f){
            double s=0;
            for (int i=1;i<=n;++i) s += stdpsi[i-1]*c[(f-1)*mdim+(i-1)];
            s=std::fabs(s);
            if (s>summax) jfill=f;
            summax=std::max(s,summax);
            if (s>0.707106781187) { ifill=f; goto done1; }
        }
        for (int f=1;f<=n;++f){
            double s=0;
            for (int i=1;i<=n;++i) s += psi[i-1]*c[(f-1)*mdim+(i-1)];
            s=std::fabs(s);
            if (s>summax) jfill=f;
            summax=std::max(s,summax);
            if (s>0.7071) { ifill=f; goto done1; }
        }
        std::printf(" CAUTION !!! SUM IN SWAP VERY SMALL, SUMMAX =%f JFILL=%d\n",summax,jfill);
        ifill=jfill;
    }
done1:
    if (ifill <= nocc) return;
    for (int i=1;i<=n;++i){
        double x=c[(nocc-1)*mdim+(i-1)];
        c[(nocc-1)*mdim+(i-1)]=c[(ifill-1)*mdim+(i-1)];
        c[(ifill-1)*mdim+(i-1)]=x;
    }
}
