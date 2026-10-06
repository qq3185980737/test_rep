// hybrid.cpp
#include "hybrid.h"
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "molkst_C.h"
#include "MOZYME_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "elemts_C.h"
#include "ijbo.h"
#include "rsp.h"
#include "local2.h"

using namespace molkst_C;
using namespace MOZYME_C;
using namespace common_arrays_C;
using namespace elemts_C;

void minloc(double* vecs, int nvec, int n) {
    auto V = [&](int i, int j)->double& { return vecs[(j-1)*nvec + (i-1)]; };
    double rot = 0.999;
    int i = 2;
    if (n != 2) {
        for (i = 2; i <= 4; ++i) {
            double sum = V(i,2)*V(i,2) + V(i,3)*V(i,3);
            if (sum > 0.1) goto L1000;
        }
L1010:
        {
            double sum = V(i,4)*V(i,4) + V(i,2)*V(i,2);
            sum = 1.0 / std::sqrt(sum);
            double alpha = V(i,4) * sum;
            double beta = V(i,2) * sum;
            for (int j = 1; j <= nvec; ++j) {
                double s = alpha*V(j,4) + beta*V(j,2);
                V(j,4) = -beta*V(j,4) + alpha*V(j,2);
                V(j,2) = s;
            }
            goto L1020;
        }
L1000:
        {
            double sum = V(i,2)*V(i,2) + V(i,3)*V(i,3);
            sum = 1.0 / std::sqrt(sum);
            double alpha = V(i,2) * sum;
            double beta = V(i,3) * sum;
            for (int j = 1; j <= nvec; ++j) {
                double s = alpha*V(j,2) + beta*V(j,3);
                V(j,3) = -beta*V(j,2) + alpha*V(j,3);
                V(j,2) = s;
            }
            goto L1010;
        }
    }
L1020:
    for (i = 2; i <= 4; ++i) {
        double sum = V(i,4)*V(i,4) + V(i,3)*V(i,3);
        if (sum > 0.1) goto L1030;
    }
    return;
L1030:
    {
        double sum = V(i,4)*V(i,4) + V(i,3)*V(i,3);
        sum = 1.0 / std::sqrt(sum);
        double alpha = V(i,4) * sum;
        double beta = V(i,3) * sum;
        for (int j = 1; j <= nvec; ++j) {
            double s = alpha*V(j,4) + beta*V(j,3);
            V(j,4) = -beta*V(j,4) + alpha*V(j,3);
            V(j,3) = s;
        }
        for (int l1 = 1; l1 <= 4; ++l1) {
            for (int l2 = l1+1; l2 <= 4; ++l2) {
                double a = rot;
                double b = std::sqrt(1.0 - a*a);
                for (int j = 1; j <= nvec; ++j) {
                    double s = a*V(j,l1) + b*V(j,l2);
                    V(j,l1) = -b*V(j,l1) + a*V(j,l2);
                    V(j,l2) = s;
                }
            }
        }
    }
}

void hybrid(double* catom) {
    static int nf_loc[17] = {0,1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    static int nl_loc[17] = {0,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    double root2 = std::sqrt(2.0);
    (void)root2;
    int loop = 0;
    for (int ii = 1; ii <= numat; ++ii) {
        int nb = nbonds[ii];
        if (nb > 15) {
            std::string el = atom_names[nat[ii]];
            std::fprintf(stdout, " For atom %d, a %s, number of bonds: %d\n",
                         ii, el.c_str(), nb);
            std::exit(1);
        }
        if (iorbs[ii] == 0) continue;
        if (iorbs[ii] == 1) {
            ++loop;
            catom[(loop-1)*morb + 0] = 1.0;
        } else {
            int k = ((nb+4)*(nb+5)) / 2;
            std::vector<double> a(191, 0.0);
            k = 10;
            for (int j = 1; j <= nb; ++j) {
                int jj = ibonds[j][ii];
                int jl = ijbo(jj, ii);
                int m = iorbs[jj];
                if (ibonds[j][ii] < ii) {
                    for (int l = 1; l <= 4; ++l)
                        a[k+l] = f[jl - m + 1 + m*l];
                } else {
                    for (int l = 1; l <= 4; ++l)
                        a[k+l] = f[jl+l];
                }
                k = k + 4 + j;
            }
            k = 4 + nb;
            for (int i = 1; i <= 4; ++i)
                for (int l = 1; l <= 4; ++l) {
                    int j = (i*(i-1))/2 + l;
                    a[j] = a[j] + j*1e-8;
                }
            for (int i = 5; i <= k; ++i)
                for (int l = 1; l <= 4; ++l) {
                    int j = (i*(i-1))/2 + l;
                    a[j] = a[j] + j*2e-4;
                }
            std::vector<double> eig(11, 0.0);
            std::vector<double> c(366, 0.0);   // k*k <= (4+15)^2 = 361, 1-based
            if (k == 1) { eig[1] = a[1]; c[1] = 1.0; }
            else rsp(a.data(), k, eig.data(), c.data());
            for (int i = 1; i <= k; ++i)
                if (c[(i-1)*k + 1] < 1e-14)
                    for (int j = 1; j <= k; ++j)
                        c[(i-1)*k+j] = -c[(i-1)*k+j];
            int j = 8 - k;
            if (j > 1) minloc(c.data(), k, j);
            local2(c.data(), k, 4, nf_loc, nl_loc, nb+1);
            for (int i = 1; i <= 4; ++i) {
                double rt2 = 0.0;
                for (int j2 = 1; j2 <= 4; ++j2)
                    rt2 += c[(i-1)*k+j2] * c[(i-1)*k+j2];
                rt2 = 1.0 / std::sqrt(rt2);
                for (int j2 = 1; j2 <= 4; ++j2)
                    catom[(i+loop-1)*morb + (j2-1)] = c[(i-1)*k+j2] * rt2;
            }
            k = std::min(4, iorbs[ii]);
            for (int i = 1; i <= k; ++i)
                for (int j2 = 5; j2 <= iorbs[ii]; ++j2)
                    catom[(i+loop-1)*morb + (j2-1)] = 0.0;
            for (int i = 5; i <= iorbs[ii]; ++i) {
                for (int j2 = 1; j2 <= iorbs[ii]; ++j2)
                    catom[(i+loop-1)*morb + (j2-1)] = 0.0;
                catom[(i+loop-1)*morb + (i-1)] = 1.0;
            }
            loop += iorbs[ii];
        }
    }
}

