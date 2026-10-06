// axis.cpp — C++ translation of MOPAC 2016 "axis.F90".
// Principal moments of inertia, molecular weight, rotational constants.

#include "axis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "to_screen_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace funcon_C;
#pragma warning(disable: 4459)

namespace {
// ddot: BLAS dot product over a 1-based run of doubles.
double ddot(int n, const double* x, int, const double* y, int) {
    double s = 0.0;
    for (int i = 1; i <= n; ++i) s += x[i] * y[i];
    return s;
}

// rsp: eigenvalues/eigenvectors of a 3x3 symmetric matrix packed as
// t(1)=M11, t(2)=M12, t(3)=M22, t(4)=M13, t(5)=M23, t(6)=M33.
// Eigenvalues sorted ascending in eig(1..3); eigenvectors as columns evec.
void rsp(const double t[7], int, double eig[4], double evec[4][4]) {
    double a[4][4];
    a[1][1] = t[1]; a[1][2] = t[2]; a[1][3] = t[4];
    a[2][1] = t[2]; a[2][2] = t[3]; a[2][3] = t[5];
    a[3][1] = t[4]; a[3][2] = t[5]; a[3][3] = t[6];
    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j) evec[i][j] = (i == j) ? 1.0 : 0.0;

    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = std::fabs(a[1][2]) + std::fabs(a[1][3]) + std::fabs(a[2][3]);
        if (off < 1.0e-30) break;
        for (int pp = 1; pp <= 2; ++pp) {
            for (int qq = pp + 1; qq <= 3; ++qq) {
                if (std::fabs(a[pp][qq]) < 1.0e-30) continue;
                double tau = (a[qq][qq] - a[pp][pp]) / (2.0 * a[pp][qq]);
                double sign = tau >= 0 ? 1.0 : -1.0;
                double theta = sign / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
                double c = 1.0 / std::sqrt(1.0 + theta * theta);
                double s = theta * c;
                double app = a[pp][pp];
                double aqq = a[qq][qq];
                a[pp][pp] = app - theta * a[pp][qq];
                a[qq][qq] = aqq + theta * a[pp][qq];
                a[pp][qq] = a[qq][pp] = 0.0;
                for (int i = 1; i <= 3; ++i) {
                    if (i != pp && i != qq) {
                        double aip = a[i][pp], aiq = a[i][qq];
                        a[i][pp] = a[pp][i] = c * aip - s * aiq;
                        a[i][qq] = a[qq][i] = s * aip + c * aiq;
                    }
                    double vip = evec[i][pp], viq = evec[i][qq];
                    evec[i][pp] = c * vip - s * viq;
                    evec[i][qq] = s * vip + c * viq;
                }
            }
        }
    }
    // gather and sort ascending
    int order[4] = {0, 1, 2, 3};
    double ev[4] = {0.0, a[1][1], a[2][2], a[3][3]};
    std::sort(order + 1, order + 4, [&](int x, int y) { return ev[x] < ev[y]; });
    double tmp[4][4];
    for (int i = 1; i <= 3; ++i) {
        eig[i] = ev[order[i]];
        for (int r = 1; r <= 3; ++r) tmp[r][i] = evec[r][order[i]];
    }
    for (int i = 1; i <= 3; ++i)
        for (int r = 1; r <= 3; ++r) evec[r][i] = tmp[r][i];
}
}

void axis(double& a, double& b, double& c, double evec[4][4]) {
    static double t[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    static double eig[4] = {0.0, 0.0, 0.0, 0.0};
    static bool first = true;
    static int icalcn = 0;

    if (icalcn != numcal) {
        icalcn = numcal;
        first = true;
    }

    double const1 = std::pow(10.0, 24) / fpc_10;
    double const2 = fpc_6 * fpc_10 * 1.0e16 / (8.0 * pi * pi * fpc_8);

    double sumwx = 0.0, sumwy = 0.0, sumwz = 0.0;
    std::vector<double> x(numat + 1, 0.0), y(numat + 1, 0.0), z(numat + 1, 0.0);

    if (mol_weight > 0) {
        for (int i = 1; i <= numat; ++i) {
            sumwx += atmass[i] * coord[0][i];
            sumwy += atmass[i] * coord[1][i];
            sumwz += atmass[i] * coord[2][i];
        }
    } else {
        mol_weight = mol_weight + (double)numat;
        for (int i = 1; i <= numat; ++i) {
            sumwx += coord[0][i];
            sumwy += coord[1][i];
            sumwz += coord[2][i];
        }
    }

    if (mol_weight > 0 && first) {
        std::printf("\n          MOLECULAR WEIGHT =%8.2f\n\n",
                    std::min(99999.99, mol_weight));
    }
    sumwx /= mol_weight;
    sumwy /= mol_weight;
    sumwz /= mol_weight;
    for (int i = 1; i <= numat; ++i) {
        x[i] = coord[0][i] - sumwx;
        y[i] = coord[1][i] - sumwy;
        z[i] = coord[2][i] - sumwz;
    }

    for (int i = 1; i <= 6; ++i) t[i] = (double)i * 1.0e-10;

    if (mol_weight > 0) {
        // weighted versions of y^2+z^2 etc.
        std::vector<double> wz2px2(numat + 1), wy2pz2(numat + 1), wx2py2(numat + 1);
        for (int i = 1; i <= numat; ++i) {
            wy2pz2[i] = y[i] * y[i] + z[i] * z[i];
            wz2px2[i] = z[i] * z[i] + x[i] * x[i];
            wx2py2[i] = x[i] * x[i] + y[i] * y[i];
        }
        t[1] += ddot(numat, &atmass[0], 1, &wy2pz2[0], 1);
        // t(2) = -sum mass*x*y
        {
            double s = 0;
            for (int i = 1; i <= numat; ++i) s -= atmass[i] * x[i] * y[i];
            t[2] += s;
        }
        t[3] += ddot(numat, &atmass[0], 1, &wz2px2[0], 1);
        {
            double s = 0;
            for (int i = 1; i <= numat; ++i) s -= atmass[i] * z[i] * x[i];
            t[4] += s;
        }
        {
            double s = 0;
            for (int i = 1; i <= numat; ++i) s -= atmass[i] * y[i] * z[i];
            t[5] += s;
        }
        t[6] += ddot(numat, &atmass[0], 1, &wx2py2[0], 1);
    } else {
        for (int i = 1; i <= numat; ++i) {
            t[1] += y[i] * y[i] + z[i] * z[i];
            t[2] -= x[i] * y[i];
            t[3] += z[i] * z[i] + x[i] * x[i];
            t[4] -= z[i] * x[i];
            t[5] -= y[i] * z[i];
            t[6] += x[i] * x[i] + y[i] * y[i];
        }
    }

    rsp(t, 3, eig, evec);

    if (mol_weight > 0 && first && keywrd.find("RC=") == std::string::npos) {
        std::printf("\n\n         ROTATIONAL CONSTANTS IN CM(-1)\n\n");
        for (int i = 1; i <= 3; ++i) {
            if (eig[i] < 3.0e-4) {
                eig[i] = 0.0;
                to_screen_C::rot[i] = 0.0;
            } else {
                to_screen_C::rot[i] = const2 / eig[i];
            }
        }
        for (int i = 1; i <= 3; ++i) to_screen_C::xyzmom[i] = eig[i] * const1;
        std::printf("         A =%14.8f   B =%14.8f   C =%14.8f\n\n",
                    to_screen_C::rot[1], to_screen_C::rot[2], to_screen_C::rot[3]);
        if (keywrd.find("RC=") == std::string::npos)
            std::printf("\n\n          PRINCIPAL MOMENTS OF INERTIA IN UNITS OF 10**(-40)*GRAM-CM**2\n\n");
        std::printf("         A =%14.4f   B =%14.4f   C =%14.4f\n\n",
                    to_screen_C::xyzmom[1], to_screen_C::xyzmom[2], to_screen_C::xyzmom[3]);
        c = to_screen_C::rot[1];
        b = to_screen_C::rot[2];
        a = to_screen_C::rot[3];
    }

    // make diagonal terms positive
    for (int i = 1; i <= 3; ++i) {
        if (evec[i][i] >= 0.0) continue;
        for (int r = 1; r <= 3; ++r) evec[r][i] = -evec[r][i];
    }
    // orient for chirality
    double sum = evec[1][1] * (evec[2][2] * evec[3][3] - evec[3][2] * evec[2][3]) +
                 evec[1][2] * (evec[2][3] * evec[3][1] - evec[2][1] * evec[3][3]) +
                 evec[1][3] * (evec[2][1] * evec[3][2] - evec[2][2] * evec[3][1]);
    if (sum < 0) {
        double mx = 1.0;
        int idx = 1;
        for (int j = 1; j <= 3; ++j) {
            if (evec[j][j] >= mx) continue;
            mx = evec[j][j];
            idx = j;
        }
        for (int r = 1; r <= 3; ++r) evec[r][idx] = -evec[r][idx];
    }
    if (keywrd.find(" NOREOR") == std::string::npos &&
        keywrd.find(" FORCETS") == std::string::npos) {
        for (int i = 1; i <= numat; ++i) {
            coord[0][i] = x[i];
            coord[1][i] = y[i];
            coord[2][i] = z[i];
        }
    }
    if (mol_weight > 0) first = false;
}
