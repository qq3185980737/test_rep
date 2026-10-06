// sort.cpp
#include "sort.h"
void sort(float* val, std::complex<float>* vec, int n) {
    for (int i=1;i<=n;++i) {
        float x=1e9f; int k=i;
        for (int j=i;j<=n;++j) if (val[j]<x) { k=j; x=val[j]; }
        for (int j=1;j<=n;++j) {
            std::complex<float> s=vec[(k-1)*n+(j-1)];
            vec[(k-1)*n+(j-1)]=vec[(i-1)*n+(j-1)];
            vec[(i-1)*n+(j-1)]=s;
        }
        val[k]=val[i];
        val[i]=x;
    }
}
