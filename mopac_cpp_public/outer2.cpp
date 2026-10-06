// outer2.cpp
#include "outer2.h"
#include <cmath>
#include "molkst_C.h"
#include "parameters_C.h"
using namespace molkst_C;
using namespace parameters_C;
extern "C" void reppd_(int* ni, int* nj, double* rij, double* ri, double* gab);
void outer2(int ni, int nj, const double* xi, const double* xj,
            double* w, int& kr, double* e1b, double* e2a,
            double& enuc, int mode, bool direct) {
    double ri[23];
    double gab = 0.0;
    if (mode == 0) {
        double x[3];
        x[0] = xi[0]-xj[0]; x[1] = xi[1]-xj[1]; x[2] = xi[2]-xj[2];
        double rij = x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
        double rijx = std::sqrt(rij);
        rij = rijx;
        reppd_(&ni, &nj, &rij, ri, &gab);
        double a = 1.0/rijx;
        x[0]*=a; x[1]*=a; x[2]*=a;
        if (std::abs(x[2]) > 0.99999999) x[2] = (x[2]>=0?1.0:-1.0);
        bool si = (natorb[ni] > 1);
        bool sj = (natorb[nj] > 1);
        double w1 = ri[0];
        int ki = 1;
        for (int i=0;i<45;++i){e1b[i]=0;e2a[i]=0;}
        e1b[0] = -w1*tore[nj];
        e2a[0] = -w1*tore[ni];
        if (!direct) w[0] = w1;
        if (sj) {
            double w2 = -ri[4]*x[0], w3 = -ri[4]*x[1], w4 = -ri[4]*x[2];
            e2a[1] = -w2*tore[ni];
            e2a[2] = -w1*tore[ni];
            e2a[3] = -w3*tore[ni];
            e2a[5] = -w1*tore[ni];
            e2a[6] = -w4*tore[ni];
            e2a[9] = -w1*tore[ni];
            if (natorb[nj] > 4) { e2a[14]=e2a[0]; e2a[20]=e2a[0]; e2a[27]=e2a[0]; e2a[35]=e2a[0]; e2a[44]=e2a[0]; }
            ki = 4;
            if (!direct) { w[1]=w2; w[2]=w3; w[3]=w4; }
        }
        if (si) {
            if (sj) {
                double w5 = -ri[1]*x[0], w6 = -ri[1]*x[1], w7 = -ri[1]*x[2];
                ki = 7;
                e1b[1] = -w5*tore[nj];
                e1b[2] = -w1*tore[nj];
                e1b[3] = -w6*tore[nj];
                e1b[5] = -w1*tore[nj];
                e1b[6] = -w7*tore[nj];
                e1b[9] = -w1*tore[nj];
                if (natorb[ni] > 4) { e1b[14]=e1b[0]; e1b[20]=e1b[0]; e1b[27]=e1b[0]; e1b[35]=e1b[0]; e1b[44]=e1b[0]; }
                if (!direct) { w[4]=w5; w[5]=w6; w[6]=w7; }
            } else {
                double w2 = -ri[1]*x[0], w3 = -ri[1]*x[1], w4 = -ri[1]*x[2];
                ki = 4;
                e1b[1] = -w2*tore[nj];
                e1b[2] = -w1*tore[nj];
                e1b[3] = -w3*tore[nj];
                e1b[5] = -w1*tore[nj];
                e1b[6] = -w4*tore[nj];
                e1b[9] = -w1*tore[nj];
                if (natorb[ni] > 4) { e1b[14]=e1b[0]; e1b[20]=e1b[0]; e1b[27]=e1b[0]; e1b[35]=e1b[0]; e1b[44]=e1b[0]; }
                if (!direct) { w[1]=w2; w[2]=w3; w[3]=w4; }
            }
        }
        enuc = tore[ni]*tore[nj]*w1;
        if (!direct) {
            if (natorb[ni]*natorb[nj] == 0) ki = 0;
            kr += ki;
        }
    } else {
        for (int i=0;i<45;++i){e1b[i]=0;e2a[i]=0;}
        enuc = 0.0;
    }
}
