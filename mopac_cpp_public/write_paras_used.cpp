// write_paras_used.cpp — C++ translation of write_paras_used (writmo.F90).
// Print the quantum-chemical parameters used for the elements present.

#include "write_paras_used.h"

#include <cstdio>
#include <string>
#include <vector>

#include "common_arrays_C.h"
#include "molkst_C.h"
#include "parameters_C.h"

using namespace common_arrays_C;
using namespace molkst_C;
using namespace parameters_C;

static const char* kSym[] = {
    "",   "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne",
    "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti",
    "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce",
    "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu",
    "Am", "Cm", "Bk", "Mi", "XX", "FM", "MD", "CB", "++", "+", "--", "- ",
    "TV"};

void write_paras_used() {
    std::vector<bool> els(99, false);
    for (int i = 1; i <= numat; ++i) {
        int e = nat[i];
        if (e < 1) e = 1;
        if (e > 98) e = 98;
        els[e] = true;
    }

    std::printf("\n\n%16sParameters used\n", "");
    std::printf("\n%5sParameter Type  Element    Parameter\n", "");

    for (int i = 1; i <= 98; ++i) {
        if (!els[i]) continue;
        std::string s = kSym[i];
#define P(tag, val) \
        if ((val) != 0.0) std::printf("%10s%7s%9s%16.6f\n", "", tag, s.c_str(), (val));
        P("USS  ", uss[i]); P("UPP  ", upp[i]); P("UDD  ", udd[i]);
        P("BETAS", betas[i]); P("BETAP", betap[i]); P("BETAD", betad[i]);
        P("ZS   ", zs[i]); P("ZP   ", zp[i]); P("ZD   ", zd[i]);
        P("ZSN  ", zsn[i]); P("ZPN  ", zpn[i]); P("ZDN  ", zdn[i]);
        P("ALP  ", alp[i]);
        P("GSS  ", gss[i]); P("GSP  ", gsp[i]); P("GPP  ", gpp[i]);
        P("GP2  ", gp2[i]); P("HSP  ", hsp[i]); P("POC  ", pocord[i]);
        if (!main_group[i] && f0sd[i] != 0.0) P("F0SD ", f0sd[i]);
        if (!main_group[i] && g2sd[i] != 0.0) P("G2SD ", g2sd[i]);
#undef P
        for (int j = 1; j <= 4; ++j) {
#define P3(tag, arr) \
            if ((arr)[i][j] != 0.0) std::printf("%12s%3s%d%10s%2s%16.6f\n", "", tag, j, "", s.c_str(), (arr)[i][j]);
            P3("FN1", guess1); P3("FN2", guess2); P3("FN3", guess3);
#undef P3
        }
        for (int j = 1; j <= i; ++j) {
            if (!els[j]) continue;
            if (alpb[i][j] > 0.01 || alpb[i][j] < -0.01)
                std::printf("%12sALPB_%s%7s%16.6f\n", "", kSym[j], kSym[i], alpb[i][j]);
            if (xfac[i][j] > 0.01 || xfac[i][j] < -0.01)
                std::printf("%12sXFAC_%s%7s%16.6f\n", "", kSym[j], kSym[i], xfac[i][j]);
        }
    }
}
