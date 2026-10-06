// bonds_for_MOZYME.cpp — C++ translation of MOPAC 2016
// "bonds_for_MOZYME.F90". Printing is stubbed; the bond-sum and valency
// logic plus the descending selection sort are translated faithfully.

#include "bonds_for_MOZYME.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "ijbo.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

namespace {
void vecprt(const std::vector<double>&, int) {}
}

void bonds_for_MOZYME() {
    bool lall = (keywrd.find(" ALLBOND") != std::string::npos);
    double sumlim = 0.01;
    if (lall) sumlim *= 0.1;

    int mtx, nper;
    if (maxtxt != 0) {
        mtx = maxtxt + 23;
        nper = 3;
    } else {
        mtx = 22;
        nper = 4;
    }
    bool mozyme_style = true;
    if (keywrd.find(" IRC") != std::string::npos ||
        keywrd.find(" DRC") != std::string::npos ||
        keywrd.find(" MINI") != std::string::npos)
        mozyme_style = false;
    if (!mozyme_style) sumlim = -1.0;  // F90: report all in classic style

    int bab_size = mozyme_style ? 100 : (nl_atoms * (nl_atoms + 1)) / 2;
    std::vector<double> bab(bab_size + 1, 0.0);
    std::vector<int> ibab(101, 0);
    std::vector<bool> l_atom_store(nl_atoms + 1, false);

    std::printf(" BOND ORDERS BETWEEN ATOMS\n");

    int newl = 0, l = 0;
    lall = true;
    for (int i = 1; i <= numat; ++i) {
        if (!mozyme_style && !l_atom[i]) continue;
        ++newl;
        double valenc = 0.0;
        int io = MOZYME_C::iorbs[i];
        int j_lim = mozyme_style ? numat : i - 1;
        if (mozyme_style) l = 0;

        if (nat[i] != 1 || lall) {
            if (io == 1) {
                int kk = ijbo(i, i) + 1;
                valenc = 2.0 * p[kk] - p[kk] * p[kk];
            } else {
                int kk = ijbo(i, i);
                for (int j = 1; j <= 4; ++j) {
                    for (int k = 1; k <= j; ++k) {
                        ++kk;
                        valenc -= p[kk] * p[kk];
                    }
                    valenc += 2.0 * p[kk];
                }
            }
            for (int j = 1; j <= j_lim; ++j) {
                if (!mozyme_style && !l_atom[j]) continue;
                if (nat[j] != 1 || lall) {
                    int jo = MOZYME_C::iorbs[j];
                    if (i != j && ijbo(i, j) >= 0) {
                        int kl = ijbo(i, j) + 1;
                        int ku = kl + io * jo - 1;
                        double sum = 0.0;
                        for (int k = kl; k <= ku; ++k) sum += p[k] * p[k];
                        if (sum > sumlim) {
                            ++l;
                            bab[l] = sum;
                            if (mozyme_style) ibab[l] = j;
                        }
                    } else if (!mozyme_style) {
                        ++l;
                        bab[l] = 0.0;
                    }
                }
            }
            if (mozyme_style) {
                if (l != 0) {
                    // Selection sort, descending bond order.
                    for (int j = 1; j <= l; ++j) {
                        double sum = 0.0;
                        int m = j;
                        for (int k = j; k <= l; ++k)
                            if (bab[k] > sum) {
                                m = k;
                                sum = bab[k] * (1.0 + 1e-4);
                            }
                        std::swap(ibab[j], ibab[m]);
                        std::swap(bab[j], bab[m]);
                    }
                    int ju = std::min(nper, l);
                    std::printf("atom %d valenc=%.3f bonds:", i, valenc);
                    for (int j = 1; j <= ju; ++j)
                        std::printf(" %s=%.3f",
                                    elemnt[nat[ibab[j]]].c_str(), bab[j]);
                    std::printf("\n");
                }
            } else {
                ++l;
                bab[l] = valenc;
            }
        }
    }
    if (!mozyme_style) {
        for (int k = 1; k <= nl_atoms; ++k) l_atom_store[k] = l_atom[k];
        for (int k = 1; k <= nl_atoms; ++k) l_atom[k] = true;
        vecprt(bab, newl);
        for (int k = 1; k <= nl_atoms; ++k) l_atom[k] = l_atom_store[k];
    }
}
