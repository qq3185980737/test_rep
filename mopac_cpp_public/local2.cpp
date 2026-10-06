// local2.cpp
#include "local2.h"
#include <cmath>
#include <vector>

void local2(double* c, int mdim, int nmos, const int* nfirst,
            const int* nlast, int numat) {
    int n = nlast[numat];
    std::vector<double> psi1(n+1), psi2(n+1);
    int niter = (nmos == 4) ? 19 : 1;
    for (int loop = 1; loop <= niter; ++loop) {
        double sum = 0.0;
        for (int i = 1; i <= nmos; ++i) {
            for (int j = 1; j <= nmos; ++j) {
                if (j == i) continue;
                double xijjj=0, xjiii=0, xiiii=0, xjjjj=0, xijij=0, xiijj=0;
                for (int k=1;k<=n;++k){ psi1[k]=c[k-1+(i-1)*mdim]; psi2[k]=c[k-1+(j-1)*mdim]; }
                for (int k1=1;k1<=numat;++k1) {
                    double dij=0,dii=0,djj=0;
                    for (int k=nfirst[k1];k<=nlast[k1];++k) {
                        dij += psi1[k]*psi2[k];
                        dii += psi1[k]*psi1[k];
                        djj += psi2[k]*psi2[k];
                    }
                    xijjj += dij*djj;
                    xjiii += dij*dii;
                    xiiii += dii*dii;
                    xjjjj += djj*djj;
                    xijij += dij*dij;
                    xiijj += dii*djj;
                }
                double aij = xijij - (xiiii+xjjjj-2.0*xiijj)/4.0;
                double bij = xjiii - xijjj;
                double ca = std::sqrt(aij*aij+bij*bij);
                double sa = aij + ca;
                if (sa > 1e-14) {
                    ca = (1.0+std::sqrt((1.0-aij/ca)/2.0))/2.0;
                    sa = std::sqrt(1.0-ca);
                    sum += std::fabs(sa);
                    ca = std::sqrt(ca);
                    for (int k=1;k<=n;++k) {
                        double p1=psi1[k],p2=psi2[k];
                        c[k-1+(i-1)*mdim] = ca*p1 + sa*p2;
                        c[k-1+(j-1)*mdim] = -sa*p1 + ca*p2;
                    }
                }
            }
        }
        if (sum < 1e-5) break;
    }
}
