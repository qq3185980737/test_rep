// rotlmo.cpp — rotate LMO p functions.
#include "rotlmo.h"
#include <vector>
#include "MOZYME_C.h"
#include "molkst_C.h"
using namespace MOZYME_C;
using molkst_C::nelecs;
using molkst_C::norbs;
using MOZYME_C::iorbs;

void rotlmo(const double* rotvec) {
    int nocc = nelecs/2;
    for (int i=1;i<=nocc;++i) {
        int loop = ncocc[i]+1;
        for (int jj=nncf[i]+1;jj<=nncf[i]+ncf[i];++jj) {
            int nj = iorbs[icocc[jj]];
            if (nj==4) {
                int ka = loop;
                double vec[4]={0,0,0,0};
                for (int j=1;j<=3;++j){
                    double sum=0;
                    for (int k=1;k<=3;++k) sum += cocc[ka+k]*rotvec[(k-1)*3+(j-1)];
                    vec[j]=sum;
                }
                for (int j=1;j<=3;++j) cocc[ka+j]=vec[j];
            }
            loop += nj;
        }
    }
    int nvir = norbs-nocc;
    for (int i=1;i<=nvir;++i) {
        int loop = ncvir[i]+1;
        for (int jj=nnce[i]+1;jj<=nnce[i]+nce[i];++jj) {
            int nj = iorbs[icvir[jj]];
            if (nj==4) {
                int ka = loop;
                double vec[4]={0,0,0,0};
                for (int j=1;j<=3;++j){
                    double sum=0;
                    for (int k=1;k<=3;++k) sum += cvir[ka+k]*rotvec[(k-1)*3+(j-1)];
                    vec[j]=sum;
                }
                for (int j=1;j<=3;++j) cvir[ka+j]=vec[j];
            }
            loop += nj;
        }
    }
}
