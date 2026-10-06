// symr.cpp — C++ translation of MOPAC 2016 "symr.F90" (complete, incl. symp).
// Build the full set of point-group symmetry operations R(9,*) and the atomic
// mapping ipo(numat,120) from the user-supplied classes in elem.
#include "symr.h"
#include "symmetry_C.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "mopend.h"
#include <cmath>
#include <cstdio>
#include <vector>

using symmetry_C::r;
using symmetry_C::ipo;
using symmetry_C::elem;
using symmetry_C::nclass;
using symmetry_C::nsym;
using symmetry_C::nent;
using symmetry_C::name;
using molkst_C::numat;
using common_arrays_C::coord;

namespace {

// symp: expand the existing operations to the full set (max 120).
void symp_impl() {
    const int maxfun = 120;
    const double tol = 1e-2;
    int i = 2, j = 1;
    while (true) {
        ++j;
        if (j > nsym) { j = 2; ++i; if (i > nsym) break; }
        if (nsym == maxfun) break;
        const int n1 = nsym + 1;
        r[1][n1] = r[1][i] * r[1][j] + r[2][i] * r[4][j] + r[3][i] * r[7][j];
        r[2][n1] = r[1][i] * r[2][j] + r[2][i] * r[5][j] + r[3][i] * r[8][j];
        r[3][n1] = r[1][i] * r[3][j] + r[2][i] * r[6][j] + r[3][i] * r[9][j];
        r[4][n1] = r[4][i] * r[1][j] + r[5][i] * r[4][j] + r[6][i] * r[7][j];
        r[5][n1] = r[4][i] * r[2][j] + r[5][i] * r[5][j] + r[6][i] * r[8][j];
        r[6][n1] = r[4][i] * r[3][j] + r[5][i] * r[6][j] + r[6][i] * r[9][j];
        r[7][n1] = r[7][i] * r[1][j] + r[8][i] * r[4][j] + r[9][i] * r[7][j];
        r[8][n1] = r[7][i] * r[2][j] + r[8][i] * r[5][j] + r[9][i] * r[8][j];
        r[9][n1] = r[7][i] * r[3][j] + r[8][i] * r[6][j] + r[9][i] * r[9][j];
        bool dup = false;
        for (int n = 1; n <= nsym; ++n) {
            double res = 0.0;
            for (int m = 1; m <= 9; ++m) res += std::fabs(r[m][n] - r[m][n1]);
            if (res < tol) { dup = true; break; }
        }
        if (dup) continue;  // not unique: discard and continue searching
        ++nsym;
        for (int n = 1; n <= numat; ++n) ipo[n][nsym] = ipo[ipo[n][j]][i];
    }
    std::printf("    FOR POINT-GROUP %s THERE ARE %3d UNIQUE SYMMETRY FUNCTIONS.\n",
                name.c_str(), nsym);
}

}  // namespace

void symr() {
    const int maxent = 6;
    ipo.assign(numat + 1, std::vector<int>(121, 0));
    bool prob = false;
    int nvalue = 0;
    // First operation is always the identity E.
    r[1][1] = 1; r[2][1] = 0; r[3][1] = 0;
    r[4][1] = 0; r[5][1] = 1; r[6][1] = 0;
    r[7][1] = 0; r[8][1] = 0; r[9][1] = 1;
    // Center the molecule.
    double x = 0, y = 0, z = 0;
    for (int i = 1; i <= numat; ++i) {
        x += coord[0][i]; y += coord[1][i]; z += coord[2][i];
        ipo[i][1] = i;
    }
    const double xa = x / numat, ya = y / numat, za = z / numat;
    for (int i = 1; i <= numat; ++i) {
        coord[0][i] -= xa; coord[1][i] -= ya; coord[2][i] -= za;
    }
    nent = 1;
    nsym = 0;
    int ndflt = 1;
    while (true) {
        ++nsym;
        ++ndflt;
        if (ndflt <= nclass) {
            // Copy symmetry operation from elem.
            int k = 0;
            nvalue = 1;
            for (int i = 1; i <= 3; ++i) {
                for (int kk = 1; kk <= 3; ++kk) r[k + kk][1 + nent] = elem[i][kk][ndflt];
                k += 3;
            }
            ++nent;
            const int n = nent;
            // Compute ipo of this function: map each atom under R.
            for (int i = 1; i <= numat; ++i) {
                const double xx = coord[0][i] * r[1][n] + coord[1][i] * r[2][n] + coord[2][i] * r[3][n];
                const double yy = coord[0][i] * r[4][n] + coord[1][i] * r[5][n] + coord[2][i] * r[6][n];
                const double zz = coord[0][i] * r[7][n] + coord[1][i] * r[8][n] + coord[2][i] * r[9][n];
                ipo[i][n] = 0;
                for (int j = 1; j <= numat; ++j) {
                    const double dist = std::fabs(xx - coord[0][j]) + std::fabs(yy - coord[1][j]) +
                                        std::fabs(zz - coord[2][j]);
                    if (dist >= 0.6) continue;
                    if (ipo[i][n] == 0) { ipo[i][n] = j; }
                    else {
                        std::printf("  ONE ATOM MAPS ONTO TWO DIFFERENT ATOMIC CENTERS\n"
                                    "  ADD KEYWORD ' NOSYM' AND RE-RUN\n");
                        prob = true;
                        break;
                    }
                }
                if (ipo[i][n] != 0) continue;
                std::printf("  ONE ATOM MAPS ONTO NO OTHER ATOM \n"
                            "  ADD KEYWORD ' NOSYM' AND RE-RUN\n");
                prob = true;
                break;
            }
        }
        if (!(nvalue != 0 && nsym < maxent)) break;
    }
    if (prob) {
        std::printf(" PROBLEM IN SYMR\n");
        mopend("PROBLEM IN SYMR");
        return;
    }
    nsym = nent;
    // Expand the existing operators to the full set.
    symp_impl();
}
