// jcarin.cpp — C++ translation of MOPAC 2016 "jcarin.F90" (complete).
// Jacobian d(Cartesian)/d(internal), by finite difference. Input: xparam
// (internal coordinates), step, preci (2-point vs 1-point). Output: b
// (iu-il+1, 3*numat*l123) Jacobian, step times too large (columns are
// Cartesian displacements per unit internal step).
#include "jcarin.h"

#include <vector>

#include "common_arrays_C.h"
#include "gmetry.h"
#include "molkst_C.h"
#include "symmetry_C.h"
#include "symtry.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace symmetry_C;

void jcarin(double* xparam, double step, bool preci, double* b, int& ncol,
            int il, int iu) {
    ncol = 3 * numat * l123;
    const int nrow = iu - il + 1;
    std::vector<double> coold(ncol + 1, 0.0);
    // Internal coordinates of the central point.
    for (int ivar = 1; ivar <= nvar; ++ivar)
        geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar];
    if (id == 0) {
        // Molecular system.
        int jvar = 0;
        for (int ivar = il; ivar <= iu; ++ivar) {
            ++jvar;
            // Step forward.
            geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar] + step;
            if (ndep != 0) symtry();
            gmetry(geo, coord);
            int j = 0;
            for (int j1 = 1; j1 <= numat; ++j1)
                for (int j2 = 1; j2 <= 3; ++j2)
                    b[(j++) * nrow + (jvar - 1)] = coord[j2 - 1][j1];
            geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar];
        }
        if (preci) {
            jvar = 0;
            for (int ivar = il; ivar <= iu; ++ivar) {
                ++jvar;
                // Step backward.
                geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar] - step;
                if (ndep != 0) symtry();
                gmetry(geo, coord);
                int j = 0;
                for (int j1 = 1; j1 <= numat; ++j1)
                    for (int j2 = 1; j2 <= 3; ++j2)
                        b[(j++) * nrow + (jvar - 1)] -= coord[j2 - 1][j1];
                geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar];
            }
        } else {
            // Central point.
            if (ndep != 0) symtry();
            gmetry(geo, coord);
            jvar = 0;
            for (int ivar = il; ivar <= iu; ++ivar) {
                ++jvar;
                int j = 0;
                for (int j1 = 1; j1 <= numat; ++j1)
                    for (int j2 = 1; j2 <= 3; ++j2)
                        b[(j++) * nrow + (jvar - 1)] -= coord[j2 - 1][j1];
            }
        }
    } else {
        // Solid state: replicate cells over the l1u/l2u/l3u ranges.
        int jvar = 0;
        for (int ivar = il; ivar <= iu; ++ivar) {
            ++jvar;
            geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar] + step;
            if (ndep != 0) symtry();
            gmetry(geo, coord);
            int ij = 0;
            for (int ii = 1; ii <= numat; ++ii)
                for (int im = -l1u; im <= l1u; ++im)
                    for (int jl = -l2u; jl <= l2u; ++jl)
                        for (int kl = -l3u; kl <= l3u; ++kl) {
                            for (int c = 1; c <= 3; ++c)
                                b[(ij + c - 1) * nrow + (jvar - 1)] =
                                    coord[c - 1][ii] + tvec[c][1] * im +
                                    tvec[c][2] * jl + tvec[c][3] * kl;
                            ij += 3;
                        }
            geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar];
        }
        if (preci) {
            jvar = 0;
            for (int ivar = il; ivar <= iu; ++ivar) {
                ++jvar;
                geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar] - step;
                if (ndep != 0) symtry();
                gmetry(geo, coord);
                int ij = 0;
                for (int ii = 1; ii <= numat; ++ii)
                    for (int im = -l1u; im <= l1u; ++im)
                        for (int jl = -l2u; jl <= l2u; ++jl)
                            for (int kl = -l3u; kl <= l3u; ++kl) {
                                for (int c = 1; c <= 3; ++c)
                                    b[(ij + c - 1) * nrow + (jvar - 1)] -=
                                        coord[c - 1][ii] + tvec[c][1] * im +
                                        tvec[c][2] * jl + tvec[c][3] * kl;
                                ij += 3;
                            }
                geo[loc[2][ivar]][loc[1][ivar]] = xparam[ivar];
            }
        } else {
            if (ndep != 0) symtry();
            gmetry(geo, coord);
            int ij = 0;
            for (int ii = 1; ii <= numat; ++ii)
                for (int im = -l1u; im <= l1u; ++im)
                    for (int jl = -l2u; jl <= l2u; ++jl)
                        for (int kl = -l3u; kl <= l3u; ++kl) {
                            coold[3 * ij] =
                                coord[0][ii] + tvec[1][1] * im +
                                tvec[1][2] * jl + tvec[1][3] * kl;
                            coold[3 * ij + 1] =
                                coord[1][ii] + tvec[2][1] * im +
                                tvec[2][2] * jl + tvec[2][3] * kl;
                            coold[3 * ij + 2] =
                                coord[2][ii] + tvec[3][1] * im +
                                tvec[3][2] * jl + tvec[3][3] * kl;
                            ij += 1;
                        }
            jvar = 0;
            for (int ivar = il; ivar <= iu; ++ivar) {
                ++jvar;
                for (int j = 0; j < ncol; ++j)
                    b[j * nrow + (jvar - 1)] -= coold[j];
            }
        }
    }
}
