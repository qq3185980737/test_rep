// symtry.cpp — C++ translation of MOPAC 2016 "symtry.F90" (complete).
// Symmetry: compute dependent bond lengths/angles/dihedrals from the
// reference internal coordinates. haddon evaluates one dependency function.
#include "symtry.h"

#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"
#include "mopend.h"
#include "symmetry_C.h"

using namespace common_arrays_C;
using namespace funcon_C;
using namespace molkst_C;
using namespace symmetry_C;

namespace {
// Evaluate one dependency function m for atom loc; returns the coordinate
// index l and value w to be stored into the dependent atom.
void haddon(double& w, int& l, int m, int loc, double fact) {
    if (m > 19 || m < 1) {
        (void)chanel_C::iw;
        mopend("UNDEFINED SYMMETRY FUNCTION USED");
        return;
    }
    const int i = loc;
    if (na[i] == 0) {
        // Cartesian relationships.
        switch (m) {
            case 1: l = 1; w = geo[i][1]; return;  // X = X
            case 2: l = 2; w = geo[i][2]; return;  // Y = Y
            case 3: l = 3; w = geo[i][3]; return;  // Z = Z
            case 4: l = 1; w = -geo[i][1]; return; // X = -X
            case 5: l = 2; w = -geo[i][2]; return; // Y = -Y
            case 6: l = 3; w = -geo[i][3]; return; // Z = -Z
            case 7: l = 1; w = geo[i][2]; return;  // X = Y
            case 8: l = 2; w = geo[i][3]; return;  // Y = Z
            case 9: l = 3; w = geo[i][1]; return;  // Z = X
            case 10: l = 1; w = -geo[i][2]; return; // X = -Y
            case 11: l = 2; w = -geo[i][3]; return; // Y = -Z
            case 12: l = 3; w = -geo[i][1]; return; // Z = -X
            case 13: l = 1; w = geo[i][3]; return;  // X = Z
            case 14: l = 2; w = geo[i][1]; return;  // Y = X
            case 15: l = 3; w = geo[i][2]; return;  // Z = Y
            case 16: l = 1; w = -geo[i][3]; return; // X = -Z
            case 17: l = 2; w = -geo[i][1]; return; // Y = -X
            case 18:
            case 19: l = 3; w = -geo[i][2]; return; // Z = -Y
            default: break;
        }
    } else {
        // Bond-angle / dihedral relationships.
        switch (m) {
            case 2: l = 2; w = geo[i][2]; return;  // set bond-angles equal
            case 3: w = geo[i][3]; break;          // set dihedrals equal
            case 4: w = pi / 2.0 - geo[i][3]; break;
            case 5: w = pi / 2.0 + geo[i][3]; break;
            case 6: w = 2.0 * pi / 3.0 - geo[i][3]; break;
            case 7: w = 2.0 * pi / 3.0 + geo[i][3]; break;
            case 8: w = pi - geo[i][3]; break;
            case 9: w = pi + geo[i][3]; break;
            case 10: w = 4.0 * pi / 3.0 - geo[i][3]; break;
            case 11: w = 4.0 * pi / 3.0 + geo[i][3]; break;
            case 12: w = 3.0 * pi / 2.0 - geo[i][3]; break;
            case 13: w = 3.0 * pi / 2.0 + geo[i][3]; break;
            case 14: w = -geo[i][3]; break;
            case 15: l = 1; w = geo[i][1] / 2.0; return;
            case 16: l = 2; w = geo[i][2] / 2.0; return;
            case 17: l = 2; w = pi - geo[i][2]; return;
            case 18:
            case 19: l = 1; w = geo[i][1] * fact; return;
            default: goto len1;
        }
        l = 3;
        return;
    }
len1:  // set bond-lengths equal
    l = 1;
    w = geo[i][1];
}
}  // namespace

void symtry() {
    int n = 0;
    for (int i = 1; i <= ndep; ++i) {
        double value = 0.0;
        int locn = 1;
        if (idepfn[i] == 19 && depmul[n + 1] > 1e-20) {
            ++n;
            haddon(value, locn, idepfn[i], locpar[i], depmul[n]);
        } else {
            haddon(value, locn, idepfn[i], locpar[i], depmul[i]);
        }
        const int j = locdep[i];
        geo[locn][j] = value;
    }
}
