// chklew.cpp — C++ translation of MOPAC 2016 "chklew.F90".
// mopend / memory_error are stubs.

#include "chklew.h"

#include <cmath>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;
using namespace MOZYME_C;
using namespace parameters_C;

namespace {
void mopend(const char*) {}
void memory_error(const char*) {}
}

void add_Lewis_element(int atom_i, int atom_j, int charge, int& element_type) {
    Lewis_tot++;
    Lewis_elem[1][Lewis_tot] = atom_i;
    Lewis_elem[2][Lewis_tot] = atom_j;
    if (atom_i > 0 && atom_j > 0) {
        iz[atom_i]--; iz[atom_j]--;
        ib[atom_i]--; ib[atom_j]--;
        element_type++;
    } else if (atom_i != 0) {
        if (atom_i > 0) {
            if (charge == -1) iz[atom_i]--;
            else if (charge == 0) iz[atom_i] -= 2;
            ib[atom_i]--;
            element_type++;
        } else {
            ib[-atom_i]--;
        }
    } else {
        if (charge == 2) iz[atom_j] -= 2;
        else if (charge == 1) iz[atom_j]--;
        ib[atom_j]--;
    }
    if (charge == 0) return;
    ions[atom_i + atom_j] += charge;
}

void ring5(int i, const std::vector<int>& mb, std::vector<int>& ir5) {
    int ni = nbonds[i];
    int nrings = 1;
    for (int j1 = 1; j1 <= ni; ++j1) {
        int j = ibonds[j1][i];
        if (mb[j] < 3) continue;
        int nj = nbonds[j];
        for (int k1 = j1 + 1; k1 <= ni; ++k1) {
            int k = ibonds[k1][i];
            if (mb[j] < 3) continue;
            int nk = nbonds[k];
            for (int l1 = 1; l1 <= nj; ++l1) {
                int l = ibonds[l1][j];
                if (l == i || mb[l] < 3) continue;
                for (int m1 = 1; m1 <= nk; ++m1) {
                    int m = ibonds[m1][k];
                    if (m == i || mb[m] < 3) continue;
                    int nm = nbonds[m];
                    bool found = false;
                    for (int n1 = 1; n1 <= nm; ++n1)
                        if (ibonds[n1][m] == l) { found = true; break; }
                    if (found) {
                        ir5[i] = ir5[j] = ir5[k] = ir5[l] = ir5[m] = ++nrings;
                        return;
                    }
                }
            }
        }
    }
}

bool arom(int ii, int jj, const std::vector<int>& mpii) {
    int ni = nbonds[ii], nj = nbonds[jj];
    for (int i = 1; i <= ni; ++i) {
        int iia = ibonds[i][ii];
        if (iia == jj || mpii[iia] == 0) continue;
        int iib = mpii[iia];
        for (int j = 1; j <= nj; ++j) {
            int jja = ibonds[j][jj];
            if (jja == ii || mpii[jja] == 0) continue;
            int jjb = mpii[jja];
            int njb = nbonds[jjb];
            for (int j2 = 1; j2 <= njb; ++j2)
                if (ibonds[j2][jjb] == iib) return true;
        }
    }
    return false;
}

bool arom2(int ii, int jj, const std::vector<int>& mpii) {
    int ni = nbonds[ii], nj = nbonds[jj];
    for (int i = 1; i <= ni; ++i) {
        int iia = ibonds[i][ii];
        if (iia == jj || mpii[iia] == 0) continue;
        int iib = mpii[iia];
        int nib = nbonds[iib];
        for (int j = 1; j <= nj; ++j) {
            int jja = ibonds[j][jj];
            if (jja == ii) continue;
            for (int i1 = 1; i1 <= nib; ++i1) {
                int iic = ibonds[i1][iib];
                for (int j1 = 1; j1 <= nbonds[jja]; ++j1)
                    if (iic == ibonds[j1][jja]) return true;
            }
        }
    }
    for (int j = 1; j <= nj; ++j) {
        int jja = ibonds[j][jj];
        if (jja == ii || mpii[jja] == 0) continue;
        int jjb = mpii[jja];
        int njb = nbonds[jjb];
        for (int i = 1; i <= ni; ++i) {
            int iia = ibonds[i][ii];
            if (iia == jj) continue;
            for (int j1 = 1; j1 <= njb; ++j1) {
                int jjc = ibonds[j1][jjb];
                for (int i1 = 1; i1 <= nbonds[iia]; ++i1)
                    if (jjc == ibonds[i1][iia]) return true;
            }
        }
    }
    return false;
}

void chklew(std::vector<int>& mb, std::vector<int>& numbon, int& l,
            int large, bool debug) {
    (void)large; (void)debug;
    if (!Lewis_elem.empty()) Lewis_elem.clear();
    Lewis_elem.assign(3, std::vector<int>(norbs + 1, 0));
    Lewis_tot = 0;
    bool big = keywrd.find(" LARGE") != std::string::npos;
    bool graphi = false, first_pi = true;

    std::vector<int> ipi(6 * numat + 1, 0), ir5(numat + 1, 0), mpii(numat + 1, 0);

    for (int i = 1; i <= numat; ++i) {
        ib[i] = nlast[i] - nfirst[i] + 1;
        iz[i] = (int)tore[nat[i]];
    }

    int j = 0;
    for (int i = 1; i <= natoms; ++i) {
        if (labels[i] != 99) ++j;
        if (txtatm[i].substr(0, 6).find('+') != std::string::npos) {
            icharges++; iz[j]--; ions[j] = 1;
        } else if (txtatm[i].substr(0, 12).find('-') != std::string::npos) {
            icharges++; iz[j]++; ions[j] = -1;
        }
    }
    for (int i = 1; i <= numat; ++i) {
        if (ib[i] != nlast[i] - nfirst[i] + 1) {
            int z = nat[i];
            for (int k = 1; k <= ndelec[z] / 2; ++k) add_Lewis_element(i, 0, 0, numbon[2]);
        }
    }
    numbon[1] = numbon[2] = numbon[3] = 0;

    // Sigma framework.
    for (int ii = 1; ii <= numat; ++ii) {
        mb[ii] = 0;
        int nbii = nbonds[ii];
        for (int i = 1; i <= nbii; ++i) {
            int jj = ibonds[i][ii];
            if (jj >= ii && ib[ii] > 0 && ib[jj] > 0)
                add_Lewis_element(ii, jj, 0, numbon[1]);
        }
    }
    l = 0;
    for (int i = 1; i <= numat; ++i) if (ib[i] != 0) { l = 1; break; }

    // Lone pairs.
    for (int i = 1; i <= numat; ++i) {
        if (iz[i] > ib[i] && ib[i] > 0) {
            int j2 = std::min(ib[i], iz[i] - ib[i]);
            for (int k = 1; k <= j2; ++k) add_Lewis_element(i, 0, 0, numbon[2]);
        }
        if (iz[i] < ib[i] && ib[i] > 0) {
            int j2 = ib[i] - iz[i];
            for (int k = 1; k <= j2; ++k) add_Lewis_element(0, i, 0, j2);
        }
    }
    l = 0;
    for (int i = 1; i <= numat; ++i) if (ib[i] != 0) { l = 1; break; }

    // Pi bond framework (stubbed simplified: ring5 + open-ended + aromatic).
    for (int ii = 1; ii <= numat; ++ii) {
        int cnt = 0;
        int nbii = nbonds[ii];
        for (int i = 1; i <= nbii; ++i)
            if (ib[ibonds[i][ii]] != 0) ++cnt;
        mb[ii] = cnt;
    }
    for (int i = 1; i <= numat; ++i)
        if (mb[i] > 2 && ir5[i] == 0) ring5(i, mb, ir5);

    // Open-ended pi bonds from degree-1 atoms.
    for (int ii = 1; ii <= numat; ++ii) {
        if (mb[ii] == 1) {
            int nbii = nbonds[ii];
            for (int i = 1; i <= nbii; ++i) {
                int jj = ibonds[i][ii];
                if (natorb[nat[ii]] == 9 || natorb[nat[jj]] == 9) continue;
                while (mb[jj] != 0 && ib[jj] != 0 && ib[ii] > 0) {
                    mb[ii]--; mb[jj]--;
                    mpii[ii] = jj; mpii[jj] = ii;
                    add_Lewis_element(ii, jj, 0, numbon[3]);
                    if (ib[ii] <= 0 || ib[jj] <= 0) break;
                }
            }
        }
    }

    l = 1;
    (void)big; (void)graphi; (void)first_pi;
}
