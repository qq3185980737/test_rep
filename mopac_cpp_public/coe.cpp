// coe.cpp
#include "coe.h"
#include <cmath>
#include <algorithm>

void coe(double x2, double y2, double z2, int norbi, int norbj,
         double* c, double& r) {
    const double rt34 = 0.86602540378444;
    const double rt13 = 0.57735026918963;
    double xy = x2*x2 + y2*y2;
    r = std::sqrt(xy + z2*z2);
    xy = std::sqrt(xy);
    double ca, cb, sa, sb;
    if (xy >= 1e-10) {
        ca = x2/xy; cb = z2/r; sa = y2/xy; sb = xy/r;
    } else {
        if (z2 <= 0.0) {
            if (z2 != 0.0) { ca=-1; cb=-1; sa=0; sb=0; goto done; }
            ca=0; cb=0; sa=0; sb=0; goto done;
        }
        ca=1; cb=1; sa=0; sb=0;
    }
done:
    for (int i=0;i<75;++i) c[i]=0.0;
    c[36] = 1.0;  // Fortran c(37) 1-based -> C++ 0-based
    int nij = std::max(norbi,norbj);
    if (nij >= 2) {
        c[55] = ca*cb;
        c[40] = ca*sb;
        c[25] = -sa;
        c[52] = -sb;
        c[37] = cb;
        c[22] = 0.0;
        c[49] = sa*cb;
        c[34] = sa*sb;
        c[19] = ca;
        if (nij >= 5) {
            double c2a = 2*ca*ca - 1.0;
            double c2b = 2*cb*cb - 1.0;
            double s2a = 2*sa*ca;
            double s2b = 2*sb*cb;
            c[74] = c2a*cb*cb + 0.5*c2a*sb*sb;
            c[59] = 0.5*c2a*s2b;
            c[44] = rt34*c2a*sb*sb;
            c[29] = -s2a*sb;
            c[14] = -s2a*cb;
            c[71] = -0.5*ca*s2b;
            c[56] = ca*c2b;
            c[41] = rt34*ca*s2b;
            c[26] = -sa*cb;
            c[11] = sa*sb;
            c[68] = rt13*sb*sb*1.5;
            c[53] = -rt34*s2b;
            c[38] = cb*cb - 0.5*sb*sb;
            c[65] = -0.5*sa*s2b;
            c[50] = sa*c2b;
            c[35] = rt34*sa*s2b;
            c[20] = ca*cb;
            c[5] = -ca*sb;
            c[62] = s2a*cb*cb + 0.5*s2a*sb*sb;
            c[47] = 0.5*s2a*s2b;
            c[32] = rt34*s2a*sb*sb;
            c[17] = c2a*sb;
            c[2] = c2a*cb;
        }
    }
}
