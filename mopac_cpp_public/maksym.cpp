// maksym.cpp — C++ translation of MOPAC 2016 "maksym.F90".
// Construct symmetry definitions automatically: identical bond lengths/angles/
// dihedrals become symmetry related. loc is (2,nvar) column-major.
#include "maksym.h"
#include <cmath>
#include <cstdio>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "symmetry_C.h"

using namespace molkst_C;
using namespace symmetry_C;
using namespace common_arrays_C;

// loc(row, col) = loc[(col-1)*2 + (row-1)]   (F90 loc(2,nvar), column-major)
namespace {
inline int& LOC(int* loc, int row, int col) { return loc[(col - 1) * 2 + (row - 1)]; }
}

void maksym(int* loc, double* xparam, double* xstore) {
    int i0 = 0;
    for (int i = 2; i <= natoms; ++i) { if (na[i] == 0) { i0 = i; break; } }
    if (i0 < natoms) {
        std::printf(" For AUTOSYM to work, geometry must be in internal coordinates.\n");
        mopend("For AUTOSYM to work, geometry must be in internal coordinates");
    }
    const double twopi = 2.0 * funcon_C::pi;
    ndep = 0;
    for (int i = 1; i <= nvar; ++i) {
        if (LOC(loc, 2, i) == 3) {
            const int j = (int)(std::copysign(0.5, xparam[i - 1]) + xparam[i - 1] / twopi);
            xparam[i - 1] = xparam[i - 1] - j * twopi;
        }
        xstore[i - 1] = xparam[i - 1];
    }
    for (int loop = 1; loop <= nvar; ++loop) {
        if (xstore[loop - 1] < -1e4) continue;
        const double xref = xstore[loop - 1];
        const int locl = LOC(loc, 2, loop);
        for (int i = loop + 1; i <= nvar; ++i) {
            if (std::abs(xref - xstore[i - 1]) >= 1e-3 || LOC(loc, 2, i) != locl) continue;
            ++ndep;
            locpar[ndep] = LOC(loc, 1, loop);
            idepfn[ndep] = locl;
            locdep[ndep] = LOC(loc, 1, i);
            xstore[i - 1] = -1e5;
        }
        // Special, common, dihedral symmetry function (m = 14).
        for (int i = loop + 1; i <= nvar; ++i) {
            if (std::abs(xref + xstore[i - 1]) >= 1e-3 || LOC(loc, 2, i) != locl) continue;
            ++ndep;
            locpar[ndep] = LOC(loc, 1, loop);
            idepfn[ndep] = 14;
            locdep[ndep] = LOC(loc, 1, i);
            xstore[i - 1] = -1e5;
        }
    }
    int j = 0;
    for (int i = 1; i <= nvar; ++i) {
        if (xstore[i - 1] <= -1e4) continue;
        ++j;
        LOC(loc, 1, j) = LOC(loc, 1, i);
        LOC(loc, 2, j) = LOC(loc, 2, i);
        xparam[j - 1] = xparam[i - 1];  // F90: xparam(j) = xparam(i)
    }
    nvar = j;
}
