// intfc.cpp
#include "intfc.h"
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "to_screen_C.h"

using namespace molkst_C;
using namespace common_arrays_C;

extern void jcarin(double*, double, bool, double*, int&, int, int);
extern "C" void jcarin_(double* xparam, double step, bool preci, double* b,
                       int* ncol, int il, int iu) {
    jcarin(xparam, step, preci, b, *ncol, il, iu);
}
extern "C" void to_screen_(const char*);

void intfc(double* fmatrx, double* xparam, double* georef,
           int* nar, int* nbr, int* ncr) {
    nvar = 0; numat = 0; l1u = 0; l2u = 0; l3u = 0;
    int j = 0;
    for (int i = 1; i <= natoms; ++i) if (nar[i] != 0) ++j;
    bool lxyz = (j == 0);
    for (int i = 1; i <= natoms; ++i) {
        na[i] = nar[i]; nb[i] = nbr[i]; nc[i] = ncr[i];
    }
    fcint.resize(4, std::vector<double>(natoms+1, 0.0));
    for (int i = 1; i <= natoms; ++i)
        for (int d = 1; d <= 3; ++d)
            geo[d][i] = georef[(i-1)*3 + d-1];
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] == 99) continue;
        ++numat;
        int ilim = lxyz ? 3 : std::min(3, i-1);
        for (j = 1; j <= ilim; ++j) {
            ++nvar;
            loc[1][nvar] = i;
            loc[2][nvar] = j;
            xparam[nvar] = geo[j][i];
        }
    }
    int n3 = 3*numat;
    double step = 1e-7;
    double stepi = 0.5 / step;
    l123 = 1;
    std::vector<double> dumy(3*numat*l123 + 1, 0.0);
    for (int i = 1; i <= nvar; ++i) {
        j = i;
        int ncol = 0;
        jcarin_(xparam, step, true, dumy.data(), &ncol, i, j);
        for (int k = 1; k <= n3; ++k) dumy[k] *= stepi;
        double sum = 0.0;
        for (j = 1; j <= n3; ++j)
            for (int k = 1; k <= n3; ++k) {
                int idx = (j >= k) ? (j*(j-1))/2 + k : (k*(k-1))/2 + j;
                sum += dumy[j] * fmatrx[idx] * dumy[k];
            }
        fcint[loc[2][i]][loc[1][i]] = sum * 1e-5;
    }
}
