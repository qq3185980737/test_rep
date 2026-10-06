// frame.cpp — C++ translation of MOPAC 2016 "frame.F90".
// Rigid translation/rotation modes added to force matrix.
#include "frame.h"

#include <algorithm>
#include <cmath>

#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

extern void axis(double&, double&, double&, double (*const evec)[4]);
extern "C" void axis_(double& a, double& b, double& c, double* rot) {
    axis(a, b, c, reinterpret_cast<double(*const)[4]>(rot));
}

void frame(std::vector<double>& fmat, int natoms_in, int mode) {
    double x, y, z;
    double rot[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    axis_(x, y, z, &rot[0][0]);
    int n3 = natoms_in * 3;
    std::vector<std::vector<double>> vib(7, std::vector<double>(n3 + 1, 0.0));
    std::vector<std::vector<double>> coord1(4, std::vector<double>(numat + 1, 0.0));
    for (int i = 1; i <= natoms_in; ++i)
        for (int j = 1; j <= 3; ++j) {
            double s = 0.0;
            for (int k = 1; k <= 3; ++k) s += coord[k-1][i] * rot[k][j];
            coord1[j][i] = s;
        }
    int j = 0;
    double wtmass = 1.0;
    if (mode == 1) {
        for (int i = 1; i <= natoms_in; ++i) {
            wtmass = std::sqrt(atmass[i]);
            ++j;
            vib[1][j] = wtmass;
            vib[5][j] = coord1[3][i] * wtmass;
            vib[6][j] = coord1[2][i] * wtmass;
            ++j;
            vib[2][j] = wtmass;
            vib[4][j] = coord1[3][i] * wtmass;
            vib[6][j] = -coord1[1][i] * wtmass;
            ++j;
            vib[3][j] = wtmass;
            vib[4][j] = -coord1[2][i] * wtmass;
            vib[5][j] = -coord1[1][i] * wtmass;
        }
    } else {
        j = 0;
        for (int i = 1; i <= natoms_in; ++i) {
            ++j;
            vib[1][j] = wtmass;
            vib[5][j] = coord1[3][i] * wtmass;
            vib[6][j] = coord1[2][i] * wtmass;
            ++j;
            vib[2][j] = wtmass;
            vib[4][j] = coord1[3][i] * wtmass;
            vib[6][j] = -coord1[1][i] * wtmass;
            ++j;
            vib[3][j] = wtmass;
            vib[4][j] = -coord1[2][i] * wtmass;
            vib[5][j] = -coord1[1][i] * wtmass;
        }
    }
    j = 1;
    for (int i = 1; i <= natoms_in; ++i) {
        for (int k = 4; k <= 6; ++k) {
            double xv = vib[k][j];
            double yv = vib[k][j+1];
            double zv = vib[k][j+2];
            vib[k][j]   = xv*rot[1][1] + yv*rot[1][2] + zv*rot[1][3];
            vib[k][j+1] = xv*rot[2][1] + yv*rot[2][2] + zv*rot[2][3];
            vib[k][j+2] = xv*rot[3][1] + yv*rot[3][2] + zv*rot[3][3];
        }
        j += 3;
    }
    double sums[7] = {0,0,0,0,0,0,0};
    for (int i = 1; i <= n3; ++i)
        for (int k = 1; k <= 6; ++k) sums[k] += vib[k][i] * vib[k][i];
    int n_trivial_local = itemp_2;
    if (n_trivial_local == 5) {
        double mn = sums[1];
        for (int i = 2; i <= 6; ++i) mn = std::min(mn, sums[i]);
        double sum = 0.01 + mn;
        for (int i = 1; i <= 6; ++i)
            if (sums[i] > sum) sums[i] = std::sqrt(1.0 / sums[i]);
    } else {
        for (int i = 1; i <= 6; ++i)
            if (sums[i] > 0) sums[i] = std::sqrt(1.0 / sums[i]);
    }
    if (id != 0) { sums[4] = sums[5] = sums[6] = 0.0; }
    double shift[7] = {0, 41000, 42000, 43000, 44000, 45000, 46000};
    for (int i = 1; i <= 6; ++i)
        for (int ii = 1; ii <= n3; ++ii) vib[i][ii] *= sums[i];
    int l = 0;
    for (int i = 1; i <= n3; ++i)
        for (int jj = 1; jj <= i; ++jj) {
            ++l;
            double s = 0.0;
            for (int k = 1; k <= 6; ++k)
                s += vib[k][i] * shift[k] * vib[k][jj];
            fmat[l] += s;
        }
}
