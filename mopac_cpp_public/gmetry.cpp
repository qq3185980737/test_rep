// gmetry.cpp — C++ translation of MOPAC 2016 "gmetry.F90".
#include "gmetry.h"
extern "C" void bangle_(double*, int, int, int, double*);
#include <cmath>
#include <vector>
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"

using namespace molkst_C;
using namespace common_arrays_C;




void renum(std::vector<std::vector<double>>& coord,
           std::vector<int>& na, std::vector<int>& nb,
           std::vector<int>& nc, int ii, int natoms) {
    int nai = na[ii];
    int nbi = nb[ii];
    double theta = 0.7853;
    int jj = 0;
    double rmin = 1e10;
    while (true) {
        for (int i = 1; i <= ii - 1; ++i) {
            if (i != nai && i != nbi) {
                double angle;
                bangle_(&coord[0][0], nai, nbi, i, &angle);
                if (angle > 1.5707963) angle = 2.0 * asin(1.0) - angle;
                if (angle >= theta) {
                    double rab = (coord[0][nbi] - coord[0][i]) * (coord[0][nbi] - coord[0][i]) +
                                 (coord[1][nbi] - coord[1][i]) * (coord[1][nbi] - coord[1][i]) +
                                 (coord[2][nbi] - coord[2][i]) * (coord[2][nbi] - coord[2][i]);
                    if (rab < rmin) { jj = i; rmin = rab; }
                }
            }
        }
        if (jj != 0) { nc[ii] = jj; return; }
        theta *= 0.5;
        if (theta < 0.0174533) theta = 0.0;
    }
}

void gmetry(std::vector<std::vector<double>>& geo,
            std::vector<std::vector<double>>& coord) {
    static int icalcn = 0, counter = 0;
    if (std::abs(step) > 1e-4) {
        double sum = 0.0;
        for (int j = 1; j <= 3; ++j)
            for (int i = 1; i <= natoms; ++i) {
                double v = geo[j][i] - geoa[j][i];
                sum += v * v;
            }
        sum = std::sqrt(sum);
        double error = (sum - step) / sum;
        for (int j = 1; j <= 3; ++j)
            for (int i = 1; i <= natoms; ++i)
                geo[j][i] = geo[j][i] - error * (geo[j][i] - geoa[j][i]);
    }
    coord[0][1] = geo[1][1];
    coord[1][1] = geo[2][1];
    coord[2][1] = geo[3][1];
    if (natoms == 1) return;
    if (na[2] == 1) {
        coord[0][2] = coord[0][1] + geo[1][2];
        coord[1][2] = coord[1][1];
        coord[2][2] = coord[2][1];
    } else {
        coord[0][2] = geo[1][2];
        coord[1][2] = geo[2][2];
        coord[2][2] = geo[3][2];
    }
    if (natoms != 2) {
        if (na[3] == 0) {
            coord[0][3] = geo[1][3];
            coord[1][3] = geo[2][3];
            coord[2][3] = geo[3][3];
        } else {
            double ccos = std::cos(geo[2][3]);
            if (na[3] == 1) coord[0][3] = coord[0][1] + geo[1][3] * ccos;
            else coord[0][3] = coord[0][2] - geo[1][3] * ccos;
            coord[1][3] = coord[1][2] + geo[1][3] * std::sin(geo[2][3]);
            coord[2][3] = coord[2][2];
        }
        for (int i = 4; i <= natoms; ++i)
            if (na[i] == 0) {
                coord[0][i] = geo[1][i];
                coord[1][i] = geo[2][i];
                coord[2][i] = geo[3][i];
            }
        for (int i = 4; i <= natoms; ++i) {
            if (na[i] != 0) {
                double cosa = std::cos(geo[2][i]);
                int mb = nb[i], mc = na[i];
                double xb = coord[0][mb] - coord[0][mc];
                double yb = coord[1][mb] - coord[1][mc];
                double zb = coord[2][mb] - coord[2][mc];
                double rbc = xb*xb + yb*yb + zb*zb;
                if (rbc < 1e-16) return;
                rbc = 1.0 / std::sqrt(rbc);
                int ma = nc[i];
                double xa = coord[0][ma] - coord[0][mc];
                double ya = coord[1][ma] - coord[1][mc];
                double za = coord[2][ma] - coord[2][mc];
                double xyb = std::sqrt(xb*xb + yb*yb);
                int k = -1;
                if (xyb <= 0.009) {
                    std::swap(xa, za); za = -xa;
                    std::swap(xb, zb); zb = -xb;
                    xyb = std::sqrt(xb*xb + yb*yb);
                    if (xyb < 0.009) return;
                    k = 1;
                }
                double costh = xb / xyb, sinth = yb / xyb;
                double xpa = xa*costh + ya*sinth;
                double ypa = ya*costh - xa*sinth;
                double sinph = zb * rbc;
                double cosph = std::sqrt(std::abs(1.0 - sinph*sinph));
                double zqa = za*cosph - xpa*sinph;
                double yza = std::sqrt(ypa*ypa + zqa*zqa);
                double coskh, sinkh;
                if (yza >= 1e-4) { coskh = ypa/yza; sinkh = zqa/yza; }
                else { coskh = 1.0; sinkh = 0.0; }
                double sina = std::sin(geo[2][i]);
                double sind = -std::sin(geo[3][i]);
                double cosd = std::cos(geo[3][i]);
                double xd = geo[1][i] * cosa;
                double yd = geo[1][i] * sina * cosd;
                double zd = geo[1][i] * sina * sind;
                double ypd = yd*coskh - zd*sinkh;
                double zpd = zd*coskh + yd*sinkh;
                double xpd = xd*cosph - zpd*sinph;
                double zqd = zpd*cosph + xd*sinph;
                double xqd = xpd*costh - ypd*sinth;
                double yqd = ypd*costh + xpd*sinth;
                if (k >= 1) {
                    double xrd = -zqd; zqd = xqd; xqd = xrd;
                }
                coord[0][i] = xqd + coord[0][mc];
                coord[1][i] = yqd + coord[1][mc];
                coord[2][i] = zqd + coord[2][mc];
            }
        }
    }
    int k = natoms;
    while (labels[k] == 107) {
        --k;
        if (k == 0) return;
    }
    ++k;
    if (icalcn != numcal) id = 0;
    if (k <= natoms) {
        int l = 0;
        for (int i = k; i <= natoms; ++i) {
            ++l;
            tvec[1][l] = coord[0][i];
            tvec[2][l] = coord[1][i];
            tvec[3][l] = coord[2][i];
        }
        id = l;
    }
    int j = 0;
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] == 99) continue;
        ++j;
        coord[0][j] = coord[0][i];
        coord[1][j] = coord[1][i];
        coord[2][j] = coord[2][i];
    }
}
