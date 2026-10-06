// mamult.cpp
#include "mamult.h"
void mamult(const double* a, const double* b, double* c, int n, double one) {
    int l = 0;
    for (int i=1;i<=n;++i) {
        int ii = ((i-1)*i)/2;
        for (int j=1;j<=i;++j) {
            int jj = ((j-1)*j)/2;
            ++l;
            double sum = 0.0;
            for (int k=1;k<=j;++k) sum += a[ii+k]*b[jj+k];
            for (int k=j+1;k<=i;++k) sum += a[ii+k]*b[((k-1)*k)/2+j];
            for (int k=i+1;k<=n;++k) {
                int kk=(k*(k-1))/2;
                sum += a[kk+i]*b[kk+j];
            }
            c[l] = sum + one*c[l];
        }
    }
}
