// solrot.cpp
#include "solrot.h"
#include "common_arrays_C.h"
#include "parameters_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include <cmath>
#include <algorithm>
using namespace common_arrays_C;
using namespace parameters_C;
using namespace funcon_C;
using namespace molkst_C;

static double trunk(double r) {
    static int icalcn=0;
    static double bound1,bound2,c,clim,cr,cr2,range;
    if (icalcn!=numcal) {
        clower = (std::min)(cutofp*2.0/3.0, 13.0);
        bound1 = clower/cutofp;
        cupper = cutofp;
        bound2 = cupper/cutofp;
        range = bound2-bound1;
        c = -0.5*bound1*bound1*cutofp/range;
        cr = 1.0 + bound1/range;
        cr2 = -1.0/(cutofp*2*range);
        clim = c + cr*cupper + cr2*cupper*cupper;
        icalcn = numcal;
    }
    if (r>clower) {
        if (r>cupper) return clim;
        return c + cr*r + cr2*r*r;
    }
    return r;
}

static void point(double& r, int ni, int nj, double* w, int& kr,
                  double* e1b, double* e2a, double& enuc) {
    r = trunk(r);
    double ee = ev*a0/r;
    double ee1 = -ee*tore[nj];
    double ee2 = -ee*tore[ni];
    int i=natorb[ni], j=natorb[nj];
    int nii=i*(i+1)/2, njj=j*(j+1)/2;
    kr = nii*njj;
    for (int k=1;k<=kr;++k) w[k]=0;
    for (int ii=1;ii<=i;++ii)
        for (int jj=1;jj<=j;++jj)
            w[((ii*(ii+1))/2-1)*njj+(jj*(jj+1))/2] = ee;
    for (int k=1;k<=nii;++k) e1b[k]=0;
    for (int k=1;k<=njj;++k) e2a[k]=0;
    for (int ii=1;ii<=i;++ii) e1b[(ii*(ii+1))/2] = ee1;
    for (int jj=1;jj<=j;++jj) e2a[(jj*(jj+1))/2] = ee2;
    enuc = -ee1*tore[ni];
}

#include "rotate.h"
extern void to_point(double, double&, double&);
static void point(double&, int, int, double*, int&, double*, double*, double&);
void solrot(int ni, int nj, const double* xi, const double* xj, double* wj,
            double* wk, int& kr, double* e1b, double* e2a, double& enuc) {
    static int icalcn=0; static double cutof2=0;
    if (icalcn!=numcal) { icalcn=numcal; cutof2=(std::min)(196.0, (2.0/3.0*cutofp)*(2.0/3.0*cutofp)); }
    double one = 1.0;
    if (std::fabs(xi[0]-xj[0])<1e-20 && std::fabs(xi[1]-xj[1])<1e-20 && std::fabs(xi[2]-xj[2])<1e-20) one=0.5;
    double wmax[2026]={0}, wsum[2026]={0}, wbits[2026]={0};
    for (int k=1;k<=45;++k){e1b[k]=0;e2a[k]=0;}
    enuc=0;
    double xdumy[3]={0,0,0};
    int kb=0;
    for (int i=-l1u;i<=l1u;++i)
      for (int j=-l2u;j<=l2u;++j)
        for (int kk=-l3u;kk<=l3u;++kk) {
            double xjuc[3];
            for (int d=1;d<=3;++d)
                xjuc[d-1] = xj[d-1] + tvec[d][1]*i + tvec[d][2]*j + tvec[d][3]*kk - xi[d-1];
            double r = xjuc[0]*xjuc[0]+xjuc[1]*xjuc[1]+xjuc[2]*xjuc[2];
            double e1bits[46]={0}, e2bits[46]={0}, enubit=0;
            if (r>cutof2) {
                r = std::sqrt(r);
                point(r, ni, nj, wbits, kb, e1bits, e2bits, enubit);
            } else {
                kb=0;
                rotate(ni, nj, xdumy, xjuc, wbits, kb, e1bits, e2bits, enubit);
            }
            for (int k=1;k<=kb;++k) wsum[k]+=wbits[k];
            if (wmax[1]<wbits[1]) for (int k=1;k<=kb;++k) wmax[k]=wbits[k];
            for (int k=1;k<=45;++k){e1b[k]+=e1bits[k];e2a[k]+=e2bits[k];}
            enuc += enubit*one;
        }
    if (one<0.9) for (int k=1;k<=kb;++k) wmax[k]=0;
    for (int k=1;k<=kb;++k){wk[k]=wmax[k];wj[k]=wsum[k];}
    kr += kb;
}


// nddo_to_point: solrot.F90 lines 97-130.  Smooth NDDO -> point-charge transition.
// wbits/e1bits/e2bits are 1-based (F90); r in Angstroms; enubit in/out.
void nddo_to_point(double* wbits, double* e1bits, double* e2bits,
                   double& enubit, double r, int ni, int nj) {
    double const_;
    if (method_pm7) {
        double dummy;
        to_point(r, dummy, const_);
    } else {
        if (r < 3.0) return;
        const_ = std::exp(-0.025 * (r - 3.0) * (r - 3.0));
    }
    double wbits_p[2026] = {0}, e1bits_p[46] = {0}, e2bits_p[46] = {0}, enubit_p = 0;
    int kb = 0;
    point(r, ni, nj, wbits_p, kb, e1bits_p, e2bits_p, enubit_p);
    for (int k = 1; k <= kb; ++k)
        wbits[k] = const_ * wbits[k] + (1.0 - const_) * wbits_p[k];
    for (int k = 1; k <= 45; ++k) {
        e1bits[k] = const_ * e1bits[k] + (1.0 - const_) * e1bits_p[k];
        e2bits[k] = const_ * e2bits[k] + (1.0 - const_) * e2bits_p[k];
    }
    enubit = const_ * enubit + (1.0 - const_) * enubit_p;
}
