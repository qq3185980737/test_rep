// mxmt.cpp
#include "mxmt.h"
void mxmt(const double* a, int nar, const double* b, int nbr, double* c,
          int ncc) {
    for (int i=1;i<=nar;++i)
        for (int j=1;j<=ncc;++j) {
            double s=0;
            for (int k=1;k<=nbr;++k)
                s+=a[(k-1)*nar+(i-1)]*b[(k-1)*ncc+(j-1)];
            c[(j-1)*nar+(i-1)]=s;
        }
}
