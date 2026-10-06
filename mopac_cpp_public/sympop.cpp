// sympop.cpp — C++ translation of MOPAC 2016 "sympop.F90".
// Find a symmetry-related atom (via ipo) that is suitable for a transition
// dipole calculation and symmetrize the Hessian block for atom i.
#include "sympop.h"
#include "symmetry_C.h"
#include "symh.h"

using symmetry_C::nsym;
using symmetry_C::ipo;

void sympop(double* h, int i, int& iskip, double* deldip) {
    for (int j = 1; j <= nsym; ++j) {
        if (ipo[i][j] >= i) continue;
        symh(h, deldip, i, j, nullptr);
        iskip = 3;  // atom ipo(i,j) is suitable for transition dipole calc'n
        return;
    }
    iskip = 0;
}
