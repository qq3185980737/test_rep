// h1elec.cpp
#include "h1elec.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include "parameters_C.h"
#include "MOZYME_C.h"
using namespace parameters_C;
using namespace MOZYME_C;
extern "C" void diat_(int* ni, int* nj, double* xjuc, double* smat);
void h1elec(int ni, int nj, const double* xi, const double* xj, double* smat) {
    double rab = (xi[0]-xj[0])*(xi[0]-xj[0]) + (xi[1]-xj[1])*(xi[1]-xj[1]) + (xi[2]-xj[2])*(xi[2]-xj[2]);
    if (rab > cutofs || (rab > 3.24 && (ni == 102 || nj == 102))) {
        for (int i = 0; i < 81; ++i) smat[i] = 0.0;
        return;
    }
    double xjuc[3] = {xj[0]-xi[0], xj[1]-xi[1], xj[2]-xi[2]};
    diat_(&ni, &nj, xjuc, smat);
    double bi[9] = {betas[ni]*0.5, betap[ni]*0.5, betap[ni]*0.5, betap[ni]*0.5,
                    betad[ni]*0.5, betad[ni]*0.5, betad[ni]*0.5, betad[ni]*0.5, betad[ni]*0.5};
    double bj[9] = {betas[nj]*0.5, betap[nj]*0.5, betap[nj]*0.5, betap[nj]*0.5,
                    betad[nj]*0.5, betad[nj]*0.5, betad[nj]*0.5, betad[nj]*0.5, betad[nj]*0.5};
    int norbi = natorb[ni], norbj = natorb[nj];
    for (int j = 1; j <= norbj; ++j)
        for (int i = 1; i <= norbi; ++i)
            smat[(j-1)*9 + (i-1)] *= (bi[i-1] + bj[j-1]);
}
