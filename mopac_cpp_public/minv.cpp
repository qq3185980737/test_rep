// minv.cpp
#include "minv.h"
#include <cmath>
#include <vector>

void minv(double* a, int n, double& d) {
    std::vector<int> l(n+1), m(n+1);
    d = 1.0;
    int nk = -n;
    for (int k=1;k<=n;++k) {
        nk += n;
        l[k] = k; m[k] = k;
        int kk = nk + k;
        double biga = a[kk-1];
        for (int j=k;j<=n;++j) {
            int iz = n*(j-1);
            for (int i=k;i<=n;++i) {
                int ij = iz + i;
                if (std::fabs(biga) >= std::fabs(a[ij-1])) continue;
                biga = a[ij-1];
                l[k] = i; m[k] = j;
            }
        }
        int j = l[k];
        if (j-k > 0) {
            int ki = k - n;
            for (int i=1;i<=n;++i) {
                ki += n;
                double hold = -a[ki-1];
                int ji = ki - k + j;
                a[ki-1] = a[ji-1];
                a[ji-1] = hold;
            }
        }
        int i = m[k];
        if (i-k > 0) {
            int jp = n*(i-1);
            for (int jj=1;jj<=n;++jj) {
                int jk = nk + jj;
                int ji = jp + jj;
                double hold = -a[jk-1];
                a[jk-1] = a[ji-1];
                a[ji-1] = hold;
            }
        }
        if (biga == 0.0) { d = 0.0; return; }
        for (int ii=1;ii<=n;++ii) {
            if (ii-k == 0) continue;
            int ik = nk + ii;
            a[ik-1] = a[ik-1]/(-biga);
        }
        for (int ii=1;ii<=n;++ii) {
            int ik = nk + ii;
            double hold = a[ik-1];
            int ij = ii - n;
            if (ii-k == 0) {}
            else {
                for (int jj=1;jj<=n;++jj) {
                    ij += n;
                    if (jj-k == 0) continue;
                    int kj = ij - ii + k;
                    a[ij-1] = hold*a[kj-1] + a[ij-1];
                }
            }
        }
        int kj = k - n;
        for (int jj=1;jj<=n;++jj) {
            kj += n;
            if (jj-k == 0) continue;
            a[kj-1] = a[kj-1]/biga;
        }
        d = std::max(-1e25, std::min(1e25, d));
        d = d*biga;
        a[kk-1] = 1.0/biga;
    }
    int k = n;
    while (true) {
        --k;
        if (k <= 0) break;
        int ii = l[k];
        if (ii-k > 0) {
            int jq = n*(k-1);
            int jr = n*(ii-1);
            for (int jj=1;jj<=n;++jj) {
                int jk = jq + jj;
                double hold = a[jk-1];
                int ji = jr + jj;
                a[jk-1] = -a[ji-1];
                a[ji-1] = hold;
            }
        }
        int jj = m[k];
        if (jj-k <= 0) continue;
        int ki = k - n;
        for (int iii=1;iii<=n;++iii) {
            ki += n;
            double hold = a[ki-1];
            int ji = ki - k + jj;
            a[ki-1] = -a[ji-1];
            a[ji-1] = hold;
        }
    }
}
