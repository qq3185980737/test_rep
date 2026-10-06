// xyzint.cpp — C++ translation of MOPAC 2016 "xyzint.F90".
// Cartesian -> internal coordinates. na/nb/nc (1-based) connectivity is
// derived from nearest-neighbor rules unless already set.
// xyz(3,numat) column-major, geo(3,numat) column-major output.
#include "xyzint.h"

#include <cmath>
#include <vector>

#include "bangle.h"
#include "dihed.h"
#include "funcon_C.h"

using namespace funcon_C;

namespace {
// Wrap a column-major double* (3,natoms) into 1-based vector-of-vectors.
std::vector<std::vector<double>> wrap_xyz(const double* xyz, int natoms) {
    std::vector<std::vector<double>> x(4, std::vector<double>(natoms + 4, 0.0));
    for (int i = 1; i <= natoms; ++i)
        for (int d = 1; d <= 3; ++d) x[d][i] = xyz[(i - 1) * 3 + (d - 1)];
    return x;
}
}  // namespace

void xyzint(double* xyz, int numat, int* na, int* nb, int* nc, double degree,
            double* geo) {
    auto X = [&](int d, int i) -> double& { return xyz[(i - 1) * 3 + (d - 1)]; };
    auto G = [&](int d, int i) -> double& { return geo[(i - 1) * 3 + (d - 1)]; };
    std::vector<std::vector<double>> xyz_v = wrap_xyz(xyz, numat);
    for (int i = 0; i < 3 * numat; ++i) geo[i] = 0.0;
    for (int i = 1; i <= numat; ++i) {
        int im1 = i - 1;
        int j = na[i];
        if (j == 0) {
            na[i] = 2;
            nb[i] = 3;
            nc[i] = 4;
            if (i == 1) continue;
            double sum = 1.0e30;
            int k = 0;
            for (j = 1; j <= im1; ++j) {
                double r = (X(1, i) - X(1, j)) * (X(1, i) - X(1, j)) +
                           (X(2, i) - X(2, j)) * (X(2, i) - X(2, j)) +
                           (X(3, i) - X(3, j)) * (X(3, i) - X(3, j));
                if (!(r < sum && na[j] != j && nb[j] != j)) continue;
                sum = r;
                k = j;
            }
            na[i] = k;
            j = k;
            if (i > 2) nb[i] = na[k];
            nc[i] = nb[k];
        }
        if (i > 3) {
            double sum;
            bangle(xyz_v, na[i], nb[i], nc[i], sum);
            if (sum < 1.0e-2 || std::abs(pi - sum) < 1.0e-2) {
                // Angle is zero or 180 degrees; search for an atom nearest 90.
                double r = 2.0;
                int l = 0;
                for (int k = 1; k <= im1; ++k) {
                    if (k == i || k == j) continue;
                    bangle(xyz_v, na[i], nb[i], k, sum);
                    if (std::abs(pi * 0.5 - sum) < r) {
                        r = std::abs(pi * 0.5 - sum);
                        l = k;
                    }
                    if (r < 0.5) break;
                }
                nc[i] = l;
            }
            double dih;
            dihed(xyz_v, i, j, nb[i], nc[i], dih);
            G(3, i) = dih;
        }
        G(3, i) *= degree;
        if (i > 2) {
            double ang;
            bangle(xyz_v, i, j, nb[i], ang);
            G(2, i) = ang;
        }
        G(2, i) *= degree;
        G(1, i) = std::sqrt((X(1, i) - X(1, j)) * (X(1, i) - X(1, j)) +
                            (X(2, i) - X(2, j)) * (X(2, i) - X(2, j)) +
                            (X(3, i) - X(3, j)) * (X(3, i) - X(3, j)));
    }
    na[1] = 0;
    nb[1] = 0;
    nc[1] = 0;
    if (numat > 1) {
        nb[2] = 0;
        nc[2] = 0;
        if (numat > 2) nc[3] = 0;
        na[2] = 1;
    }
}

// Vector-based convenience overload used by getpdb.cpp and friends.
void xyzint(const std::vector<double>& xyz, int natoms,
            const std::vector<int>& na, const std::vector<int>& nb,
            const std::vector<int>& nc, double degree,
            std::vector<double>& coord) {
    std::vector<double> xyz_f = xyz;
    std::vector<int> na_f = na, nb_f = nb, nc_f = nc;
    std::vector<double> coord_f = coord;
    xyzint(xyz_f.data(), natoms, na_f.data(), nb_f.data(), nc_f.data(), degree,
           coord_f.data());
    coord = coord_f;
}
