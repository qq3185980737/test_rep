// symt.cpp
#include "symt.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "mat33.h"
#include <cmath>
using symmetry_C::nsym;
using symmetry_C::ipo;
using symmetry_C::r;
using molkst_C::numat;

void symt(double* h, double* deldip, double* ha) {
    if (nsym < 2) return;
    int m = numat*(9*numat+3)/2;
    for (int i=1;i<=m;++i) ha[i]=0;
    for (int n=1;n<=nsym;++n) {
        for (int i=1;i<=numat;++i) {
            int i3=3*i, i6=6*i;
            for (int j=1;j<i;++j) {
                int k=ipo[i][n], l=ipo[j][n];
                int j3=3*j, k3=3*k, k6=6*k, l3=3*l, l6=6*l;
                double temp[10];
                int iel;
                if (k>l) {
                    iel=(k3*(k3-1))/2+l3;
                    temp[9]=h[iel]; temp[8]=h[iel-1]; temp[7]=h[iel-2];
                    temp[6]=h[iel-k3+1]; temp[5]=h[iel-k3]; temp[4]=h[iel-k3-1];
                    temp[3]=h[iel-k6+3]; temp[2]=h[iel-k6+2]; temp[1]=h[iel-k6+1];
                } else {
                    iel=(l3*(l3-1))/2+k3;
                    temp[9]=h[iel]; temp[6]=h[iel-1]; temp[3]=h[iel-2];
                    temp[8]=h[iel-l3+1]; temp[5]=h[iel-l3]; temp[2]=h[iel-l3-1];
                    temp[7]=h[iel-l6+3]; temp[4]=h[iel-l6+2]; temp[1]=h[iel-l6+1];
                }
                double temp2[10];
                double ra[10];
        for (int qq = 1; qq <= 9; ++qq) ra[qq] = r[qq][n];
        mat33(ra, temp, temp2);
                iel=(i3*(i3-1))/2+j3;
                ha[iel]=temp2[9]+ha[iel];
                ha[iel-1]=temp2[8]+ha[iel-1];
                ha[iel-2]=temp2[7]+ha[iel-2];
                ha[iel-i3+1]=temp2[6]+ha[iel-i3+1];
                ha[iel-i3]=temp2[5]+ha[iel-i3];
                ha[iel-i3-1]=temp2[4]+ha[iel-i3-1];
                ha[iel-i6+3]=temp2[3]+ha[iel-i6+3];
                ha[iel-i6+2]=temp2[2]+ha[iel-i6+2];
                ha[iel-i6+1]=temp2[1]+ha[iel-i6+1];
            }
            int k=ipo[i][n], k3=3*k, k6=6*k;
            int iel=(k3*(k3+1))/2;
            double temp[10];
            temp[9]=h[iel]; temp[8]=h[iel-1]; temp[7]=h[iel-2];
            temp[6]=temp[8]; temp[5]=h[iel-k3]; temp[4]=h[iel-k3-1];
            temp[3]=temp[7]; temp[2]=temp[4]; temp[1]=h[iel-k6+1];
            double temp2[10];
            double ra[10];
    for (int qq = 1; qq <= 9; ++qq) ra[qq] = r[qq][n];
    mat33(ra, temp, temp2);
            iel=(i3*(i3+1))/2;
            ha[iel]=temp2[9]+ha[iel];
            ha[iel-1]=temp2[8]+ha[iel-1];
            ha[iel-2]=temp2[7]+ha[iel-2];
            ha[iel-i3]=temp2[5]+ha[iel-i3];
            ha[iel-i3-1]=temp2[4]+ha[iel-i3-1];
            ha[iel-i6+1]=temp2[1]+ha[iel-i6+1];
            // dipole: F90 deldip(3, 3*numat) column-major, 1-based
            //   deldip(i,j) = deldip[(j-1)*3 + (i-1)]
            temp[9]=deldip[(k*3-1)*3+2]; temp[8]=deldip[(k*3-1)*3+1]; temp[7]=deldip[(k*3-1)*3+0];
            temp[6]=deldip[(k*3-2)*3+2]; temp[5]=deldip[(k*3-2)*3+1]; temp[4]=deldip[(k*3-2)*3+0];
            temp[3]=deldip[(k*3-3)*3+2]; temp[2]=deldip[(k*3-3)*3+1]; temp[1]=deldip[(k*3-3)*3+0];
            mat33(&r[1][n], temp, temp2);
            int di3=i*3;
            deldip[(di3-1)*3+2]+=temp2[9]; deldip[(di3-1)*3+1]+=temp2[8]; deldip[(di3-1)*3+0]+=temp2[7];
            deldip[(di3-2)*3+2]+=temp2[6]; deldip[(di3-2)*3+1]+=temp2[5]; deldip[(di3-2)*3+0]+=temp2[4];
            deldip[(di3-3)*3+2]+=temp2[3]; deldip[(di3-3)*3+1]+=temp2[2]; deldip[(di3-3)*3+0]+=temp2[1];
        }
    }
    for (int i=1;i<=m;++i) h[i]=ha[i]/nsym;
    const int nd3 = 3*numat*3;
    for (int i=0;i<nd3;++i) deldip[i]/=nsym;
}
