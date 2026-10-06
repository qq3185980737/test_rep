// dimens.cpp — C++ translation of MOPAC 2016 "dimens.F90".
// Computes the three molecular dimensions (longest axis, second-longest,
// thickness) and prints the MOLECULAR DIMENSIONS table.  Like the Fortran
// original it temporarily rotates the coordinates and restores them on exit.

#include "dimens.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;   // nat
using namespace elemts_C;          // elemnt
using namespace molkst_C;          // numat

namespace {

// symopr: coord' = R^T . coord  (F90 symopr with jump >= 0).
void symopr(int numat, std::vector<std::vector<double>>& coord, int jump,
            const double r[3][3]) {
    if (jump >= 0) {
        for (int i = 1; i <= numat; ++i) {
            double help[3] = { coord[0][i], coord[1][i], coord[2][i] };
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int jj = 0; jj < 3; ++jj) s += r[jj][j] * help[jj];
                coord[j][i] = s;
            }
        }
    } else {
        for (int i = 1; i <= numat; ++i) {
            double help[3] = { coord[0][i], coord[1][i], coord[2][i] };
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int jj = 0; jj < 3; ++jj) s += r[j][jj] * help[jj];
                coord[j][i] = s;
            }
        }
    }
}

}  // namespace

void dimens(std::vector<std::vector<double>>& coord, int iw) {
    (void)iw;
    if (numat == 1) return;

    std::vector<std::vector<double>> store(4,
        std::vector<double>(numat + 1, 0.0));
    for (int j = 0; j < 3; ++j)
        for (int i = 1; i <= numat; ++i) store[j][i] = coord[j][i];

    int kmax = 1, lmax = 1;
    double rmax = 0.0, r = 0.0;

    if (numat > 20) {
        // Pick any atom (atom 1) and walk out to the most distant atom,
        // then to the midpoint, repeating until the distance converges.
        int l = 1;
        double x1 = coord[0][l], y1 = coord[1][l], z1 = coord[2][l];
        double rabmax = 0.0;
        int loop = 0;
        for (;;) {
            ++loop;
            rmax = 0.0;
            int k = 1;
            for (int j = 1; j <= numat; ++j) {
                double rr = (x1 - coord[0][j]) * (x1 - coord[0][j]) +
                            (y1 - coord[1][j]) * (y1 - coord[1][j]) +
                            (z1 - coord[2][j]) * (z1 - coord[2][j]);
                if (rr > rmax) { rmax = rr; k = j; }
            }
            if (std::fabs(rmax - rabmax) < 1.0e-5) break;
            rabmax = std::max(rmax, rabmax);
            kmax = k;
            lmax = l;
            x1 = 0.5 * (x1 + coord[0][k]);
            y1 = 0.5 * (y1 + coord[1][k]);
            z1 = 0.5 * (z1 + coord[2][k]);
            rmax = 0.0;
            l = 1;
            for (int j = 1; j <= numat; ++j) {
                double rr = (x1 - coord[0][j]) * (x1 - coord[0][j]) +
                            (y1 - coord[1][j]) * (y1 - coord[1][j]) +
                            (z1 - coord[2][j]) * (z1 - coord[2][j]);
                if (rr > rmax) { rmax = rr; l = j; }
            }
            x1 = coord[0][l];
            y1 = coord[1][l];
            z1 = coord[2][l];
            if (loop >= 10) break;
        }
    } else {
        // Search for the absolute largest distance.
        rmax = 0.0;
        for (int k = 1; k <= numat; ++k) {
            double x1 = coord[0][k], y1 = coord[1][k], z1 = coord[2][k];
            for (int l = 1; l < k; ++l) {
                double rr = (x1 - coord[0][l]) * (x1 - coord[0][l]) +
                            (y1 - coord[1][l]) * (y1 - coord[1][l]) +
                            (z1 - coord[2][l]) * (z1 - coord[2][l]);
                if (rr > rmax) { rmax = rr; kmax = k; lmax = l; }
            }
        }
    }

    int k = kmax, l = lmax;
    double x1 = coord[0][k] - coord[0][l];
    double y1 = coord[1][k] - coord[1][l];
    double z1 = coord[2][k] - coord[2][l];
    double xy = x1 * x1 + y1 * y1;
    r = std::sqrt(xy + z1 * z1);
    xy = std::sqrt(xy);

    double ca, cb, sa, sb;
    if (xy < 1.0e-10) {
        if (z1 < 0.0)      { ca = -1.0; cb = -1.0; sa = 0.0; sb = 0.0; }
        else if (z1 > 0.0) { ca =  1.0; cb =  1.0; sa = 0.0; sb = 0.0; }
        else               { ca =  0.0; cb =  0.0; sa = 0.0; sb = 0.0; }
    } else {
        ca = x1 / xy;
        cb = z1 / r;
        sa = y1 / xy;
        sb = xy / r;
    }
    double c[3][3];
    c[0][2] = ca * cb;  c[0][1] = -sa;      c[0][0] = ca * sb;
    c[1][2] = sa * cb;  c[1][1] = ca;       c[1][0] = sa * sb;
    c[2][2] = -sb;      c[2][1] = 0.0;      c[2][0] = cb;
    symopr(numat, coord, 1, c);

    double dim[3];
    int ij[3][2];
    dim[0] = r;
    ij[0][0] = k;
    ij[0][1] = l;

    // Second dimension: most distant pair in the rotated Y-Z plane.
    int kk = k, ll = l;
    double rabmax = 0.0;
    int loop = 0;
    double y1p = coord[1][l], z1p = coord[2][l];
    r = 0.0;
    for (;;) {
        ++loop;
        if (loop > 10) { k = kk; l = ll; break; }
        rmax = 0.0;
        int jk = 1;
        for (int j = 1; j <= numat; ++j) {
            double rr = (y1p - coord[1][j]) * (y1p - coord[1][j]) +
                        (z1p - coord[2][j]) * (z1p - coord[2][j]);
            if (rr > rmax) { rmax = rr; jk = j; }
        }
        k = jk;
        if (std::fabs(rmax - rabmax) < 1.0e-5) {
            r = std::sqrt(rmax + 1.0e-20);
            break;
        }
        if (rmax > rabmax) {
            rabmax = rmax;
            kk = k;
            ll = l;
        }
        y1p = 0.5 * (y1p + coord[1][k]);
        z1p = 0.5 * (z1p + coord[2][k]);
        rmax = 0.0;
        int jl = 1;
        for (int j = 1; j <= numat; ++j) {
            double rr = (y1p - coord[1][j]) * (y1p - coord[1][j]) +
                        (z1p - coord[2][j]) * (z1p - coord[2][j]);
            if (rr > rmax) { rmax = rr; jl = j; }
        }
        l = jl;
        y1p = coord[1][l];
        z1p = coord[2][l];
    }
    double y1v = coord[1][k] - coord[1][l];
    double z1v = coord[2][k] - coord[2][l];
    r = std::sqrt(y1v * y1v + z1v * z1v + 1.0e-20);
    ca = y1v / r;
    sa = z1v / r;
    for (int i = 1; i <= numat; ++i)
        coord[2][i] = -sa * coord[1][i] + ca * coord[2][i];

    if (r > dim[0] + 1.0e-10) {
        // This should not normally happen; keep dim(1) and assign dim(2).
        dim[1] = r;
        ij[1][0] = k;
        ij[1][1] = l;
    } else {
        dim[1] = r;
        ij[1][0] = k;
        ij[1][1] = l;
    }

    // Third dimension: extent along the new Z axis.
    double ymin = 1.0e16, ymax = -1.0e16;
    int k3 = 1, l3 = 1;
    for (int i = 1; i <= numat; ++i) {
        if (coord[2][i] > ymax) { k3 = i; ymax = coord[2][i]; }
        if (coord[2][i] < ymin) { l3 = i; ymin = coord[2][i]; }
    }
    dim[2] = ymax - ymin;
    ij[2][0] = k3;
    if (k3 == l3) l3 = (l3 == 1) ? 2 : 1;
    ij[2][1] = l3;

    std::printf("\n%9sMOLECULAR DIMENSIONS (Angstroms)\n\n", "");
    std::printf("%9s   Atom       Atom       Distance\n", "");
    for (int j = 0; j < 3; ++j)
        std::printf("%11s%-2s%6d%3s%-2s%6d%12.5f\n", "",
                    elemts_C::elemnt[nat[ij[j][0]]].c_str(), ij[j][0],
                    "", elemts_C::elemnt[nat[ij[j][1]]].c_str(), ij[j][1],
                    dim[j]);

    // Restore coordinates.
    for (int j = 0; j < 3; ++j)
        for (int i = 1; i <= numat; ++i) coord[j][i] = store[j][i];
}
