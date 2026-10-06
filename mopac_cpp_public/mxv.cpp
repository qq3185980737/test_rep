// mxv.cpp
#include "mxv.h"
void mxv(const double* a, int nar, const double* vecx, int nbr, double* vecy) {
    for (int i=1;i<=nar;++i) {
        double s=0;
        for (int k=1;k<=nbr;++k)
            s+=a[(k-1)*nar+(i-1)]*vecx[k];
        vecy[i]=s;
    }
}
