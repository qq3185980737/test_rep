// dhcore.cpp — C++ translation of MOPAC 2016 "dhcore.F90".
// 1- and 2-electron integral derivatives w.r.t. coord(natx,nati) by
// 2-point finite difference. h/ww 1-based packed; e1b/e2a 0-based (rotate
// output layout); di/ddi 0-based column-major 9x9.
#include "dhcore.h"
#include <cmath>
#include <vector>
#include "common_arrays_C.h"
#include "molkst_C.h"
#include "h1elec.h"
#include "rotate.h"
using namespace common_arrays_C;
using namespace molkst_C;

void dhcore(double* coord, double* h, double* ww, double& enuclr,
            int nati, int natx, double step) {
    int N = molkst_C::norbs * (molkst_C::norbs + 1) / 2;
    for (int i = 1; i <= N; ++i) h[i] = 0.0;   // F90 h(1:N), 1-based-padding (h[0] unused)
    enuclr = 0.0;
    int kr = 1;
    int ia = nfirst[nati], ib = nlast[nati], ni = nat[nati];
    // work copy of coordinates (input stays const)
    std::vector<double> c(coord, coord + 3 * molkst_C::numat);
    double csave = c[(nati - 1) * 3 + (natx - 1)];
    for (int j = 1; j <= molkst_C::numat; ++j) {
        if (j == nati) continue;
        int ja = nfirst[j], jb = nlast[j], nj = nat[j];
        double di[81] = {}, ddi[81] = {};
        c[(nati - 1) * 3 + (natx - 1)] = csave + step;
        h1elec(ni, nj, &c[(nati - 1) * 3], &c[(j - 1) * 3], di);
        c[(nati - 1) * 3 + (natx - 1)] = csave - step;
        h1elec(ni, nj, &c[(nati - 1) * 3], &c[(j - 1) * 3], ddi);
        // fill atom-other-atom one-electron block
        int i2 = 0;
        if (ia > ja) {
            for (int i1 = ia; i1 <= ib; ++i1) {
                int ij = (i1 * (i1 - 1)) / 2 + ja - 1;
                i2++;
                for (int m = 0; m < jb - ja + 1; ++m)
                    h[ij + 1 + m] += di[(i2 - 1) + m * 9] - ddi[(i2 - 1) + m * 9];
            }
        } else {
            for (int i1 = ja; i1 <= jb; ++i1) {
                int ij = (i1 * (i1 - 1)) / 2 + ia - 1;
                i2++;
                // F90: h(ij+1:...) += di(1:ib-ia+1, i2) - ddi(1:ib-ia+1, i2)
                // di is column-major: di(row, col) physical = (col-1)*9 + (row-1)
                // row = O AO (m+1), col = H AO (i2)  =>  physical = (i2-1)*9 + m
                for (int m = 0; m < ib - ia + 1; ++m)
                    h[ij + 1 + m] += di[(i2 - 1) * 9 + m] - ddi[(i2 - 1) * 9 + m];
            }
        }
        int kro = kr;
        double wjd[2027] = {}, dwjd[2027] = {};
        double e1b[45] = {}, de1b[45] = {}, e2a[45] = {}, de2a[45] = {};
        double enuc = 0.0, denuc = 0.0;
        c[(nati - 1) * 3 + (natx - 1)] = csave + step;
        rotate(ni, nj, &c[(nati - 1) * 3], &c[(j - 1) * 3], wjd, kr, e1b, e2a, enuc);
        kr = kro;
        c[(nati - 1) * 3 + (natx - 1)] = csave - step;
        rotate(ni, nj, &c[(nati - 1) * 3], &c[(j - 1) * 3], dwjd, kr, de1b, de2a, denuc);
        int k = kr - kro;
        if (nati == 1 && j == 2) { fprintf(stderr, "[DH] k=%d kro=%d kr=%d\n", k, kro, kr); fflush(stderr); }
        if (kr > 0) {
            for (int i = 0; i <= k; ++i) wjd[i] -= dwjd[i];   // F90 wjd(1:k+1): 0-based physical
            for (int j = 0; j < k; ++j) ww[kro + j] = wjd[j];  // F90 do i=0,k-1: ww(i+kro)=wjd(i+1)
            if (nati == 1 && natx == 1 && j == 2) {
                fprintf(stderr, "[DH] OH1 wjd0..9:"); for (int q = 0; q < 10; ++q) fprintf(stderr, " %+.4e", wjd[q]); fprintf(stderr, "\n"); fflush(stderr);
                fprintf(stderr, "[DH] OH1 dwj0..9:"); for (int q = 0; q < 10; ++q) fprintf(stderr, " %+.4e", dwjd[q]); fprintf(stderr, "\n"); fflush(stderr);
                fprintf(stderr, "[DH] OH1 ww1..10:"); for (int q = 1; q <= 10; ++q) fprintf(stderr, " %+.4e", ww[kro - 1 + q]); fprintf(stderr, "\n"); fflush(stderr);
            }
            if (nati == 1 && natx == 1 && j == 2) {
                fprintf(stderr, "[DH] OH1 e1bd0..9:"); for (int q = 0; q < 10; ++q) fprintf(stderr, " %+.4e", e1b[q] - de1b[q]); fprintf(stderr, "\n"); fflush(stderr);
                fprintf(stderr, "[DH] OH1 e2ad0..1:"); for (int q = 0; q < 2; ++q) fprintf(stderr, " %+.4e", e2a[q] - de2a[q]); fprintf(stderr, "\n"); fflush(stderr);
                fprintf(stderr, "[DH] OH1 h1d0..3:"); for (int q = 0; q < 4; ++q) fprintf(stderr, " %+.4e", di[q * 9 + 0] - ddi[q * 9 + 0]); fprintf(stderr, "\n"); fflush(stderr);
            }
        }
        c[(nati - 1) * 3 + (natx - 1)] = csave;
        enuclr += enuc - denuc;
        // electron-nuclear attraction: atom i (e1b) and atom j (e2a) blocks
        i2 = 0;
        for (int i1 = ia; i1 <= ib; ++i1) {
            int ii = (i1 * (i1 - 1)) / 2 + ia - 1;
            int cnt = i1 - ia + 1;
            for (int m = 0; m < cnt; ++m)
                h[ii + 1 + m] += e1b[i2 + m] - de1b[i2 + m];
            i2 += cnt;
        }
        i2 = 0;
        for (int i1 = ja; i1 <= jb; ++i1) {
            int ii = (i1 * (i1 - 1)) / 2 + ja - 1;
            int cnt = i1 - ja + 1;
            for (int m = 0; m < cnt; ++m)
                h[ii + 1 + m] += e2a[i2 + m] - de2a[i2 + m];
            i2 += cnt;
        }
    }
}
