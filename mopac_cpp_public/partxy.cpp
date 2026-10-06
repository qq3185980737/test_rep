// partxy.cpp — C++ translation of MOPAC 2016 "partxy.F90".
// PQ34(PQ) = <P,Q|C3,C4> (MNDO two-electron transform).
// c34/pq34/w are 1-based arrays (element 0 padding); formxy from M12.
#include "partxy.h"
#include <vector>
#include "common_arrays_C.h"
#include "formxy.h"
#include "molkst_C.h"

using common_arrays_C::nfirst;
using common_arrays_C::nlast;
using molkst_C::numat;

namespace {
// Fortran nb(-1:8) = {0,1,0,0,10,0,0,0,0,45}; stored with +1 offset.
const int nb[10] = {0, 1, 0, 0, 10, 0, 0, 0, 0, 45};
}  // namespace

void partxy(const double* c34, double* pq34, const double* w) {
    for (int i = 0; i < molkst_C::lm61; ++i) pq34[i] = 0.0;
    if (molkst_C::numat == 2 && nfirst[1] == 1 && nlast[1] == 1) {
        fprintf(stderr, "[PXY] c34[0..2]=%+.4f %+.4f %+.4f\n", c34[0], c34[1], c34[2]); fflush(stderr);
    }
    int iband = 1, kr = 1, ls = 0;
    for (int ii = 1; ii <= numat; ++ii) {
        int ia = nfirst[ii], ib = nlast[ii];
        if (ib >= ia) {
            ls += iband;
            iband = nb[ib - ia + 1];  // nb(ib-ia), Fortran offset -1 -> +1
            int lsp = ls - 1;
            int jband = 1, js = 0;
            for (int jj = 1; jj <= ii - 1; ++jj) {
                js += jband;
                jband = nb[nlast[jj] - nfirst[jj] + 1];
                if (jband != 0) {
                    // formxy (M12) uses 1-based vector views + kr advance.
                    // w is 1-based-padding: element kr lives at physical kr, so
                    // the view must start at w+kr (NOT w+kr-1, which is element kr-1).
                    std::vector<double> wv(w + kr, w + kr + iband * jband + 1);
                    std::vector<double> ca(c34 + ls - 1, c34 + ls - 1 + iband + 1);
                    std::vector<double> cb(c34 + js - 1, c34 + js - 1 + jband + 1);
                    // formxy outputs are 1-based vectors; align slot t to 1-based pq34(ls-1+t).
                    std::vector<double> wca(iband + 1, 0.0);
                    std::vector<double> wcb(jband + 1, 0.0);
                    for (int t = 1; t <= iband; ++t) wca[t] = pq34[ls - 1 + t - 1];
                    for (int t = 1; t <= jband; ++t) wcb[t] = pq34[js - 1 + t - 1];
                    int kr2 = kr;
                    formxy(wv, kr2, wca, wcb, ca, cb, iband, jband);
                    // formxy outputs are 1-based vectors (element 0 padding).
                    for (int t = 1; t <= iband; ++t) pq34[ls - 1 + (t - 1)] = wca[t];
                    for (int t = 1; t <= jband; ++t) pq34[js - 1 + (t - 1)] = wcb[t];
                    kr = kr2;
                }
            }
            lsp = ls - 1;
            for (int i = ia; i <= ib; ++i) {
                double aa = 1.0;
                for (int j = ia; j <= i; ++j) {
                    if (i == j) aa = 0.5;
                    lsp += 1;
                    double sum = 0.0;
                    int lsw = ls - 1;
                    for (int k = ia; k <= ib; ++k) {
                        double bb = 1.0;
                        for (int l = ia; l <= k; ++l) {
                            if (l == k) bb = 0.5;
                            lsw += 1;
                            sum += c34[lsw] * w[kr] * bb;
                            kr += 1;
                        }
                    }
                    pq34[lsp - 1] += sum * aa;  // lsp is 1-based (Fortran)
                }
            }
        }
    }
    if (molkst_C::numat == 2 && nfirst[1] == 1 && nlast[1] == 1) {
        fprintf(stderr, "[PXY] pq34[0..2]=%+.4f %+.4f %+.4f\n", pq34[0], pq34[1], pq34[2]); fflush(stderr);
    }
}
