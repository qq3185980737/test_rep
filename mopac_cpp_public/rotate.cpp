// rotate.cpp — C++ translation of MOPAC 2016 "rotate.F90".
// F90 semantics: w is 1-based (w(1) = first integral). rotatd_ writes 0-based,
// so results are shifted to w[1..kr-1] here.
#include "rotate.h"
#include <cmath>
#include <cstring>
#include "molkst_C.h"
#include "parameters_C.h"
using namespace molkst_C;
using namespace parameters_C;
extern "C" void rotatd_(int* ni, int* nj, const double* xi, const double* xj, double* w, int* kr, double* enuc);
extern "C" void elenuc_(int*, int*, int*, int*, double* en);
extern "C" void nddo_to_point_(double* w, double* e1b, double* e2a, double* enuc, double* rij, int* ni, int* nj);
void rotate(int ni, int nj, const double* xi, const double* xj,
            double* w, int& kr, double* e1b, double* e2a, double& enuc) {
    double x[3];
    x[0] = xi[0]-xj[0]; x[1] = xi[1]-xj[1]; x[2] = xi[2]-xj[2];
    double rij = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    if (rij < 0.00002) {
        for (int i=0;i<45;++i){e1b[i]=0;e2a[i]=0;}
        for (int i=1;i<2026;++i) w[i]=0;
        enuc = 0.0;
        return;
    }
    double rijx = std::sqrt(rij);
    rij = rijx;
    int li = natorb[ni], lj = natorb[nj];
    double w0[2026] = {};
    int k0 = kr;              // hcore passes the absolute kr position; keep it
    int kr2 = 1;              // rotatd_ advances from 1 (relative segment length)
    rotatd_(&ni, &nj, xi, xj, w0, &kr2, &enuc);
    for (int k = 1; k < kr2; ++k) w[k-1] = w0[k-1];   // w[0..] = global w[k0..k0+len-1]
    kr = k0 + (kr2 - 1);      // absolute advance by segment length

    double en[171] = {};
    int ilow = 1, iup = li, jlow = li + 1, jup = li + lj;   // F90: elenuc(1, li, li+1, li+lj, en)
    elenuc_(&ilow, &iup, &jlow, &jup, en);
    int ik = 0;
    for (int i=1;i<=li;++i)
      for (int k=1;k<=i;++k) {
        ik++;
        e1b[ik-1] = en[(i*(i-1))/2 + k - 1];
      }
    ik = 0;
    for (int i=li+1;i<=li+lj;++i)
      for (int k=li+1;k<=i;++k) {
        ik++;
        e2a[ik-1] = en[(i*(i-1))/2 + k - 1];
      }
    if (id == 3 && !method_pm7) {
        double w1[2026] = {};
        for (int k = 0; k < kr - k0; ++k) w1[k] = w[k];
        nddo_to_point_(w1, e1b, e2a, &enuc, &rij, &ni, &nj);
        for (int k = 0; k < kr - k0; ++k) w[k] = w1[k];
    }
}
