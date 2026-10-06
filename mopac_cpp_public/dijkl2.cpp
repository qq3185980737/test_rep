// dijkl2.cpp
#include "dijkl2.h"
#include "meci_C.h"
#include "molkst_C.h"
using meci_C::nmos;
using meci_C::dijkl;
using meci_C::xy;
using molkst_C::norbs;
namespace {
double ddot(int n,const double* a,const double* b){
    double s=0; for(int i=0;i<n;++i) s+=a[i]*b[i]; return s;
}
}
void dijkl2(const std::vector<std::vector<double>>& dc) {
    int ij=0;
    for (int i=1;i<=nmos;++i)
    for (int j=1;j<=i;++j) {
        ++ij; bool lij=(i==j); int kl=0;
        for (int k=1;k<=i;++k) {
            int ll=(k==i)?j:k;
            for (int l=1;l<=ll;++l) {
                ++kl; bool lkl=(k==l);
                const double* djk = &dijkl[(kl-1)*norbs*nmos+(j-1)*norbs];
                const double* dik = &dijkl[(kl-1)*norbs*nmos+(i-1)*norbs];
                const double* dij = &dijkl[(ij-1)*norbs*nmos+(l-1)*norbs];
                const double* dkj = &dijkl[(ij-1)*norbs*nmos+(k-1)*norbs];
                double val = ddot(norbs, dc[i-1].data(), djk);
                if (lij && lkl && j==k) val*=4.0;
                else {
                    if (lij) val*=2.0;
                    else val += ddot(norbs, dc[j-1].data(), dik);
                    double val2 = ddot(norbs, dc[k-1].data(), dij);
                    if (lkl) val += val2*2.0;
                    else val += val2 + ddot(norbs, dc[l-1].data(), dkj);
                }
                auto idx=[&](int a,int b,int c,int d){
                    return ((d-1)*nmos*nmos*nmos+(c-1)*nmos*nmos+(b-1)*nmos+(a-1));
                };
                xy[idx(i,j,k,l)]=val; xy[idx(i,j,l,k)]=val;
                xy[idx(j,i,k,l)]=val; xy[idx(j,i,l,k)]=val;
                xy[idx(k,l,i,j)]=val; xy[idx(k,l,j,i)]=val;
                xy[idx(l,k,i,j)]=val; xy[idx(l,k,j,i)]=val;
            }
        }
    }
}
