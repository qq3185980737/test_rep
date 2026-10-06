// symh.cpp
#include "symh.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "mat33.h"
using symmetry_C::r;
using symmetry_C::ipo;
using molkst_C::numat;

void symh(double* h, double* dip, int i, int n, const int* ipo_arr) {
    (void)ipo_arr;
    int k = ipo[i][n];
    int i3=3*i, i6=6*i, k3=3*k, k6=6*k;
    for (int j=numat;j>=i+1;--j) {
        int l = ipo[j][n];
        int l3=3*l, l6=6*l, j3=3*j, j6=6*j;
        double temp[10], temp2[10];
        int iel;
        if (k>l) {
            iel=(k3*(k3-1))/2+l3;
            temp[9]=0.5*h[iel]; temp[8]=0.5*h[iel-1]; temp[7]=0.5*h[iel-2];
            temp[6]=0.5*h[iel-k3+1]; temp[5]=0.5*h[iel-k3]; temp[4]=0.5*h[iel-k3-1];
            temp[3]=0.5*h[iel-k6+3]; temp[2]=0.5*h[iel-k6+2]; temp[1]=0.5*h[iel-k6+1];
        } else {
            iel=(l3*(l3-1))/2+k3;
            double fact = (l<i)?0.5:1.0;
            temp[9]=fact*h[iel]; temp[6]=fact*h[iel-1]; temp[3]=fact*h[iel-2];
            temp[8]=fact*h[iel-l3+1]; temp[5]=fact*h[iel-l3]; temp[2]=fact*h[iel-l3-1];
            temp[7]=fact*h[iel-l6+3]; temp[4]=fact*h[iel-l6+2]; temp[1]=fact*h[iel-l6+1];
        }
        double ra[10];
        for (int qq = 1; qq <= 9; ++qq) ra[qq] = r[qq][n];
        mat33(ra, temp, temp2);
        iel=(j3*(j3-1))/2+i3;
        h[iel]=temp2[9]; h[iel-j3+1]=temp2[8]; h[iel-j6+3]=temp2[7];
        h[iel-1]=temp2[6]; h[iel-j3]=temp2[5]; h[iel-j6+2]=temp2[4];
        h[iel-2]=temp2[3]; h[iel-j3-1]=temp2[2]; h[iel-j6+1]=temp2[1];
    }
    int iel=(k3*(k3+1))/2;
    double temp[10], temp2[10];
    temp[9]=0.5*h[iel]; temp[8]=0.5*h[iel-1]; temp[7]=0.5*h[iel-2];
    temp[6]=temp[8]; temp[5]=0.5*h[iel-k3]; temp[4]=0.5*h[iel-k3-1];
    temp[3]=temp[7]; temp[2]=temp[4]; temp[1]=0.5*h[iel-k6+1];
    double ra[10];
    for (int qq = 1; qq <= 9; ++qq) ra[qq] = r[qq][n];
    mat33(ra, temp, temp2);
    iel=(i3*(i3+1))/2;
    h[iel]=temp2[9]; h[iel-1]=temp2[8]; h[iel-2]=temp2[7];
    h[iel-i3]=temp2[5]; h[iel-i3-1]=temp2[4]; h[iel-i6+1]=temp2[1];
    // dipole: dip(3,*) column-major: dip(i,j)=dip[(j-1)*3+(i-1)]
    temp[9]=dip[(k*3-1)*3+2]; temp[8]=dip[(k*3-1)*3+1]; temp[7]=dip[(k*3-1)*3+0];
    temp[6]=dip[(k*3-2)*3+2]; temp[5]=dip[(k*3-2)*3+1]; temp[4]=dip[(k*3-2)*3+0];
    temp[3]=dip[(k*3-3)*3+2]; temp[2]=dip[(k*3-3)*3+1]; temp[1]=dip[(k*3-3)*3+0];
    double ra2[10];
    for (int qq = 1; qq <= 9; ++qq) ra2[qq] = r[qq][n];
    mat33(ra2, temp, temp2);
    dip[(i*3-1)*3+2]=temp2[9]; dip[(i*3-1)*3+1]=temp2[8]; dip[(i*3-1)*3+0]=temp2[7];
    dip[(i*3-2)*3+2]=temp2[6]; dip[(i*3-2)*3+1]=temp2[5]; dip[(i*3-2)*3+0]=temp2[4];
    dip[(i*3-3)*3+2]=temp2[3]; dip[(i*3-3)*3+1]=temp2[2]; dip[(i*3-3)*3+0]=temp2[1];
    int im1t3=(i-1)*3;
    int istart=(im1t3*(im1t3+1))/2+1;
    for (int q=istart;q<=iel;++q) h[q]+=h[q];
}
