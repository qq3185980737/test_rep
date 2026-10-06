// xyzcry.cpp
#include "xyzcry.h"
#include <cmath>
#include <cstdio>
// tvec(3,3) column-major: tvec(i,j) = tvec[(j-1)*3 + i-1]
// dxyz(3,numat) column-major: dxyz(i,j) = dxyz[(j-1)*3 + i-1]
void xyzcry(double* tvec, int numat, double* dxyz, int iw) {
    auto tv = [&](int i, int j)->double& { return tvec[(j-1)*3 + i-1]; };
    auto dx = [&](int i, int j)->double& { return dxyz[(j-1)*3 + i-1]; };
    double sum = std::sqrt(tv(2,1)*tv(2,1) + tv(3,1)*tv(3,1));
    if (sum > 1e-6) {
        double ca = tv(2,1)/sum, sa = tv(3,1)/sum;
        for (int i=1;i<=3;++i) {
            double s1 = tv(2,i)*ca + tv(3,i)*sa;
            tv(3,i) = -tv(2,i)*sa + tv(3,i)*ca;
            tv(2,i) = s1;
        }
        for (int i=1;i<=numat;++i) {
            double s1 = dx(2,i)*ca + dx(3,i)*sa;
            dx(3,i) = -dx(2,i)*sa + dx(3,i)*ca;
            dx(2,i) = s1;
        }
        sum = std::sqrt(tv(1,1)*tv(1,1) + tv(2,1)*tv(2,1));
        ca = tv(1,1)/sum; sa = tv(2,1)/sum;
        for (int i=1;i<=3;++i) {
            double s1 = tv(1,i)*ca + tv(2,i)*sa;
            tv(2,i) = -tv(1,i)*sa + tv(2,i)*ca;
            tv(1,i) = s1;
        }
        for (int i=1;i<=numat;++i) {
            double s1 = dx(1,i)*ca + dx(2,i)*sa;
            dx(2,i) = -dx(1,i)*sa + dx(2,i)*ca;
            dx(1,i) = s1;
        }
    }
    sum = std::sqrt(tv(2,2)*tv(2,2) + tv(3,2)*tv(3,2));
    if (sum > 1e-6) {
        double ca = tv(2,2)/sum, sa = tv(3,2)/sum;
        for (int i=2;i<=3;++i) {
            double s1 = tv(2,i)*ca + tv(3,i)*sa;
            tv(3,i) = -tv(2,i)*sa + tv(3,i)*ca;
            tv(2,i) = s1;
        }
        for (int i=1;i<=numat;++i) {
            double s1 = dx(2,i)*ca + dx(3,i)*sa;
            dx(3,i) = -dx(2,i)*sa + dx(3,i)*ca;
            dx(2,i) = s1;
        }
    }
    for (int i=1;i<=3;++i) {
        double s = 0.0;
        for (int j=1;j<=i;++j) s += tv(j,i)*tv(j,i);
        for (int j=1;j<=i;++j) tv(j,i) /= s;
    }
    for (int i=1;i<=numat;++i) {
        dx(3,i) /= tv(3,3);
        dx(2,i) -= dx(3,i)*tv(2,3);
        dx(1,i) -= dx(3,i)*tv(1,3);
        dx(2,i) /= tv(2,2);
        dx(1,i) -= dx(2,i)*tv(1,2);
        dx(1,i) /= tv(1,1);
    }
    if (iw > 0) {  // F90: write(iw, '(A)') ' Fractional Unit Cell Derivatives'
        std::fprintf(stderr, " Fractional Unit Cell Derivatives\n");
        for (int i = 1; i <= numat; ++i)
            std::fprintf(stderr, "%4d%12.5f%12.5f%12.5f\n", i, dx(1,i), dx(2,i), dx(3,i));
    }
}
