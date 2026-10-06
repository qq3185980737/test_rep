// chkion.cpp — C++ translation of MOPAC 2016 "chkion.F90".
// add_Lewis_element / mopend are external stubs.

#include "chkion.h"

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

namespace {
void add_Lewis_element(int, int, int, int&) {}
void mopend(const char*) {}
int sgn(double x) { return (x > 0) ? 1 : ((x < 0) ? -1 : 0); }

// Default oxidation states (1..107).
const int ox_ref[130] = {
    0,
    1, 0,                                   // 1 H, 2 He
    1, 2, 3, 4, 3, -2, -1, 0,              // 3 Li .. 10 Ne
    1, 2, 3, 4, 3, -2, -1, 0,              // 11 Na .. 18 Ar
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, 3, -2, -1, 0,  // 19..36
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, 3, -2, -1, 0,  // 37..54
    1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,  // 55..71 (14*3)
    4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3, -2, -1, 0,           // 72..86
    1, 2, 2, 4, 2, 2, 4, 3, 2, 3, 2, 3, 2, 0, 0,             // 87..101
    0, 0, 0, 0, 0, 0, 0                                       // 102..107
};
const int ox_ref_sp[130] = {
    0,
    1, 0,
    1, 2, 3, 4, 3, -2, -1, 0,
    1, 2, 3, 4, 3, -2, -1, 0,
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, 3, -2, 1, 0,
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, 3, -2, 1, 0,
    1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2,                // 55..70
    4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3, -2, 0, 0,
    1, 2, 2, 4, 2, 2, 4, 3, 2, 3, 2, 3, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0
};
const int ox_ref_spd1[130] = {
    0,
    1, 0,
    1, 2, 3, 4, 5, 6, 3, 0,
    1, 2, 3, 4, 5, 6, 3, 0,
    1, 2, 2, 4, 5, 6, 7, 2, 3, 2, 1, 2, 3, 4, 5, 6, 3, 0,
    1, 2, 2, 4, 5, 6, 7, 2, 3, 2, 1, 2, 3, 4, 5, 6, 3, 0,
    1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2,
    4, 5, 6, 7, 2, 3, 2, 1, 2, 3, 2, 3, 4, 0, 0,
    1, 2, 2, 4, 2, 2, 4, 5, 6, 7, 2, 3, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0
};
const int ox_ref_spd2[130] = {
    0,
    1, 0,
    1, 2, 3, 4, -3, 2, -1, 0,
    1, 2, 3, 4, -5, 6, 7, 0,
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, -5, 6, 5, 0,
    1, 2, 2, 4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, -5, 6, 7, 0,
    1, 2, 2, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2,
    4, 3, 2, 3, 2, 3, 2, 1, 2, 3, 4, -3, 2, 0, 0,
    1, 2, 2, 4, 2, 2, 4, 3, 2, 3, 2, 3, 2, 0, 0,
    0, 0, 0, 0, 0, 0, 0
};
}

void chkion(std::vector<int>& ox_calc, int& n_lone_pairs,
            const std::vector<char>& atom_charge) {
    for (int i = 1; i <= numat; ++i) {
        int k = ions[i];
        int jj = (nbonds[i] == 0) ? sgn(ox_ref[nat[i]]) : sgn(ox_ref_sp[nat[i]]);
        for (int j = 1; j <= Lewis_tot; ++j) {
            if (Lewis_elem[1][j] != 0 && Lewis_elem[2][j] != 0) {
                if (Lewis_elem[1][j] == i) k += jj;
                if (Lewis_elem[2][j] == i) k += jj;
            }
        }
        ox_calc[i] = k;
        if (ib[i] == 0) {
            ions[i] = iz[i];
            iz[i] = 0;
        }
    }

    int dummy = 0;
    for (int outer_loop = 1; outer_loop <= 2; ++outer_loop) {
        for (int i = 1; i <= numat; ++i) {
            if (ib[i] <= 0) continue;
            if (nat[i] == 6) {
                if (outer_loop == 1) {
                    if (iz[i] > 1) {
                        add_Lewis_element(i, 0, 0, n_lone_pairs);
                    } else {
                        for (int m = 1; m <= nbonds[i]; ++m) {
                            int mm = ibonds[m][i];
                            int g = (int)parameters_C::tore[nat[mm]];
                            if (g == 5) { add_Lewis_element(0, i, 1, dummy); break; }
                            if (g == 6) { add_Lewis_element(i, 0, -1, n_lone_pairs); break; }
                        }
                    }
                } else {
                    add_Lewis_element(-i, 0, 0, dummy);
                }
            } else {
                int ni = nat[i];
                int jj;
                do {
                    if (atom_charge[i] == ' ') {
                        if (std::abs(ox_calc[i]) > std::abs(ox_ref_spd1[ni]))
                            jj = ox_ref_spd2[ni] - ox_calc[i];
                        else if (std::abs(ox_calc[i]) > std::abs(ox_ref_sp[ni]))
                            jj = ox_ref_spd1[ni] - ox_calc[i];
                        else if (nbonds[i] > 0)
                            jj = ox_ref_sp[ni] - ox_calc[i];
                        else
                            jj = ox_ref[ni] - ox_calc[i];
                    } else if (atom_charge[i] == '+') jj = 1;
                    else if (atom_charge[i] == '-') jj = -1;
                    else jj = 0;

                    if (jj == 0 || ib[i] == 0) break;
                    if (jj > 1) { add_Lewis_element(0, i, 2, dummy); ox_calc[i] += 2; }
                    else if (jj == 1) { add_Lewis_element(0, i, 1, dummy); ox_calc[i] += 1; }
                    else if (jj == -1) { add_Lewis_element(i, 0, -1, n_lone_pairs); ox_calc[i] -= 1; }
                    else { add_Lewis_element(i, 0, -2, n_lone_pairs); ox_calc[i] -= 2; }
                } while (true);
            }
        }
    }

    int nocc_local = 0;
    for (int i = 1; i <= Lewis_tot; ++i)
        if (Lewis_elem[1][i] != 0) ++nocc_local;
    int mm = nocc_local - nelecs / 2;
    for (int i = 1; i <= mm; ++i)
        for (int j = 1; j <= Lewis_tot; ++j) {
            int jj = Lewis_elem[1][j];
            if (jj < 0) {
                jj = -jj;
                Lewis_elem[2][j] = jj;
                Lewis_elem[1][j] = 0;
                iz[jj]--;
                for (; ib[jj] > 0;) add_Lewis_element(0, jj, 0, dummy);
                break;
            }
        }
    for (int j = 1; j <= Lewis_tot; ++j) {
        int jj = Lewis_elem[1][j];
        if (jj < 0) {
            jj = -jj;
            Lewis_elem[1][j] = jj;
            n_lone_pairs++;
            if (iz[jj] == 1) iz[jj] = 0;
            if (ib[jj] == 1) {
                add_Lewis_element(0, jj, 0, dummy);
                ib[jj] = 0;
            }
        }
    }
}
