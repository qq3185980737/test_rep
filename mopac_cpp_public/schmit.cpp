// schmit.cpp
#include "schmit.h"
#include <cmath>
void schmit(double* u, int n, int ndim) {
    const double small=0.01, one=1.0;
    int ii=0;
    for (int k=1;k<=n;++k) {
        int k1=k-1;
        double dot=0.0;
        for (int i=1;i<=n;++i) { double v=u[(i-1)+(k-1)*ndim]; dot+=v*v; }
        if (std::fabs(dot)<1e-20) goto l100;
        { double scale=one/std::sqrt(dot);
          for (int i=1;i<=n;++i) u[(i-1)+(k-1)*ndim]*=scale; }
      l30:
        if (k1==0) continue;
        int npass=0;
      l40:
        ++npass;
        for (int j=1;j<=k1;++j) {
            double d=0.0;
            for (int i=1;i<=n;++i) d+=u[(i-1)+(j-1)*ndim]*u[(i-1)+(k-1)*ndim];
            for (int i=1;i<=n;++i) u[(i-1)+(k-1)*ndim]-=d*u[(i-1)+(j-1)*ndim];
        }
        dot=0.0;
        for (int i=1;i<=n;++i) { double v=u[(i-1)+(k-1)*ndim]; dot+=v*v; }
        if (std::fabs(dot)<1e-20) goto l100;
        if (dot<small && npass>2) goto l100;
        { double scale=one/std::sqrt(dot);
          for (int i=1;i<=n;++i) u[(i-1)+(k-1)*ndim]*=scale; }
        if (dot<small) goto l40;
        continue;
      l100:
        ++ii;
        u[(ii-1)+(k-1)*ndim]=one;
        goto l30;
    }
}
