// empiri.cpp — C++ translation of "empiri" (inside writmo.F90).
// Count elements (C, H, N, O first, then others in encounter order) and
// build the empirical formula string stored in molkst_C::formula.

#include "empiri.h"

#include <string>
#include <vector>
#include <cstdio>

#include "common_arrays_C.h"
#include "molkst_C.h"

using namespace common_arrays_C;
using namespace molkst_C;

static const char* kElem[108] = {
    "",   "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne",
    "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti",
    "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce",
    "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu",
    "Am", "Cm", "Bk", "Mi", "XX", "+3", "-3", "Cb", "++", "+ ", "--", "- ",
    "Tv"};

void empiri() {
    std::vector<int> llab, mlab;
    // C, H, N, O first (atomic numbers 6,1,7,8).
    llab.push_back(6);
    llab.push_back(1);
    llab.push_back(7);
    llab.push_back(8);
    mlab.push_back(0);
    mlab.push_back(0);
    mlab.push_back(0);
    mlab.push_back(0);
    int nlab = 4;

    for (int i = 1; i <= numat; ++i) {
        int j = nat[i];
        int k = 0;
        for (int t = 1; t <= nlab; ++t)
            if (j == llab[t - 1]) { k = t; break; }
        if (k == 0) {
            ++nlab;
            llab.push_back(j);
            mlab.push_back(1);
        } else {
            ++mlab[k - 1];
        }
    }

    // Compact out zero-count entries.
    int out = 0;
    for (int i = 0; i < nlab; ++i) {
        if (mlab[i] != 0) {
            mlab[out] = mlab[i];
            llab[out] = llab[i];
            ++out;
        }
    }
    nlab = out;

    std::string f = "           Empirical Formula: ";  // F90: (10X,A,1X, ...)
    for (int i = 0; i < nlab; ++i) {
        const char* sym = kElem[llab[i]];
        // A1 for one-letter symbols ("C " -> "C"), A2 for two-letter ("Cu")
        f += (sym[1] == ' ') ? std::string(1, sym[0]) : std::string(sym);
        // count: F90 prints mlab>1 with Iw; mlab==1 is zeroed and I1.0 prints
        // a blank field (no digit) -- verified against real MOPAC output
        // "C4 H22 N10 O2 Cu Cl6".
        if (mlab[i] > 1) f += std::to_string(mlab[i]);
        f += " ";  // 1X between elements
    }
    // F90: write(formula(len_trim+1:),"(a,i6,a)") "  =",numat," atoms"
    char tail[48];
    std::snprintf(tail, sizeof(tail), " =%6d atoms", numat);
    f += tail;
    formula = f;
}
