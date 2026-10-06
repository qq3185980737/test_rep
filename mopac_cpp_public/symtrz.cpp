// symtrz.cpp — C++ translation of MOPAC 2016 "symtrz.F90".
// Determine point group and symmetrize orbitals / vibrational modes.
// makopr is a permanent gap in the 2016 tree (no source) and is supplied
// as an external stub by the caller (or real implementation when available).
#include "symtrz.h"
#include <algorithm>
#include <string>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "meci_C.h"
#include "molkst_C.h"
#include "mult33.h"
#include "symmetry_C.h"
#include "molsymy.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace meci_C;
using namespace symmetry_C;

extern void symoir(int itype, double* vects, double* eigs, int nvecs, double* r, int imat);
extern void makopr(int numat, double* cotim, int& ierror, double* r);  // permanent gap (stub by caller)

void symtrz(double* vects, double* eigs, int itype, int geteig) {
    // cotim(3,numat) column-major, copied from coord (direct indexing).
    std::vector<double> cotim(3 * numat + 1, 0.0);
    for (int i = 1; i <= numat; ++i)
        for (int j = 1; j <= 3; ++j) cotim[(i - 1) * 3 + (j - 1)] = coord[j-1][i];

    std::vector<int> nat_store(numat + 1);
    for (int i = 1; i <= numat; ++i) nat_store[i] = nat[i];

    // Fortran: if (size(jelem) < numat*20) return;   (symtrz.F90:69)
    // jelem(20,n) is allocated in setup_mopac_arrays; size(jelem)=20*n, so the
    // test is false whenever jelem is allocated for n >= numat.  C++ jelem is
    // vector<vector<int>> (21 rows x n+1 cols), so emulate the allocation guard:
    if (jelem.empty() || (int)jelem[0].size() < numat + 1) return;
    for (int i = 1; i <= numat; ++i)
        nat[i] = nat[i] * 20 + (int)(atmass[i] * 20.0 + 0.5);  // F90 Nint(atmass*20)

    double r[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    int ierror = 0;
    molsym(&cotim[0], ierror, r);
    for (int i = 1; i <= numat; ++i) nat[i] = nat_store[i];
    if (moperr) return;

    int i = std::max({norbs, maxci + 20, 3 * numat, lab});
    namo.assign(i + 1, " ");
    jndex.assign(i + 1, 0);
    for (int k = 1; k <= i; ++k) jndex[k] = k;

    if ((geteig || itype == 2) && ierror == 0) makopr(numat, &cotim[0], ierror, r);  // permanent gap: makopr

    if (itype == 2 && !geteig) {
        if (keywrd.find("SYMTRZ") != std::string::npos) {
            std::printf(" Symmetry Operations in SYMTRZ\n");
            for (int k = 1; k <= nclass; ++k) {
                std::printf(" Operation: %d\n", k);
                for (int a = 1; a <= 3; ++a)
                    std::printf("%12.6f%12.6f%12.6f\n", elem[a][1][k], elem[a][2][k], elem[a][3][k]);
            }
            std::printf(" Orientation Matrix\n");
            for (int a = 0; a < 3; ++a)
                std::printf("%12.6f%12.6f%12.6f\n", r[a * 3 + 0], r[a * 3 + 1], r[a * 3 + 2]);
        }
        for (int k = 2; k <= nclass; ++k) mult33(r, k);
    }

    if (ierror == 0 && geteig) {
        int nvecs;
        if (itype == 2) {
            for (int k = 1; k <= numat; ++k) jndex[k] = 3;
            nvecs = 3 * numat;
        } else {
            for (int k = 1; k <= numat; ++k) jndex[k] = nlast[k] - nfirst[k] + 1;
            nvecs = norbs;
        }
        if (nvecs >= 1) {
            int imat = nvecs;
            if (itype == 3) imat = lab;
            symoir(itype, vects, eigs, nvecs, r, imat);
        }
    }
}
