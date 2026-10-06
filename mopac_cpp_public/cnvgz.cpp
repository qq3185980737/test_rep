// cnvgz.cpp
#include "cnvgz.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include "molkst_C.h"
#include "MOZYME_C.h"
using namespace molkst_C;
using namespace MOZYME_C;
void cnvgz(std::vector<double>& pnew, std::vector<double>& p,
           std::vector<double>& p1, std::vector<double>& p2,
           std::vector<double>& p3, int niter,
           const std::vector<int>& idiag) {
    for (int i = 1; i <= norbs; ++i) {
        int j = idiag[i];
        p3[i] = pnew[j];
        p2[i] = p[j];
    }
    pmax = 0.0;
    for (int i = 1; i <= mpack; ++i) {
        double sa = std::abs(pnew[i] - p[i]);
        pmax = std::max(pmax, sa);
    }
    if (use_three_point_extrap) {
        if (niter % 3 == 0) {
            double faca = 0.0, facb = 0.0;
            for (int i = 1; i <= norbs; ++i) {
                double sa = std::abs(p3[i] - p2[i]);
                faca += sa * sa;
                facb += (p3[i] - 2.0 * p2[i] + p1[i]) * (p3[i] - 2.0 * p2[i] + p1[i]);
            }
            if (facb > 0.0 && faca < 100.0 * facb) {
                double fac = std::sqrt(faca / facb);
                for (int i = 1; i <= mpack; ++i) pnew[i] += fac * (pnew[i] - p[i]);
            }
        }
        if (niter > 3) {
            double damp = 0.05;
            if (pmax > damp) {
                for (int i = 1; i <= norbs; ++i) {
                    int j = idiag[i];
                    if (std::abs(p3[i] - p2[i]) > damp) {
                        double d = p3[i] - p2[i];
                        pnew[j] = p2[i] + (d >= 0 ? damp : -damp);
                        pnew[j] = std::min(2.0, std::max(pnew[j], 0.0));
                    }
                }
            }
        }
    }
    p1 = p2;
    p = pnew;
}
