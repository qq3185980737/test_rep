// mpcsyb.cpp — SYBYL output driver + Mulliken populations (F90 mpcsyb.F90).
// Output unit iw maps to stdout (project convention).
#define _CRT_SECURE_NO_WARNINGS
#include "mpcsyb.h"
#include "molkst_C.h"
#include "common_arrays_C.h"
#include "chanel_C.h"
#include "parameters_C.h"
#include "elemts_C.h"
#include "to_screen.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
using namespace molkst_C;
using namespace common_arrays_C;
using namespace chanel_C;
using namespace parameters_C;
using namespace elemts_C;

// F90 mpcpop(icok): total Mulliken populations on atoms from the diagonal
// elements of the density matrix pb.
//   icok =  1 : SYBYL present, do analysis + write SYBYL output
//   icok =  0 : SYBYL present but MULLIK absent - no analysis
//   icok = -1 : MULLIK present but SYBYL absent - analysis only
void mpcpop(int icok, FILE* isyb) {
    int i, if_ = 0, il = 0, j, k;
    std::vector<double> pop(numat + 1, 0.0);
    if (icok > -1) std::fprintf(isyb, "%4d     MULLIKEN POPULATION AND CHARGE\n", icok);
    std::vector<double> chrg(numat + 1, 0.0);   // F90: common chrg, local here
    if (icok != 0) {
        for (i = 1; i <= numat; ++i) {
            if_ = nfirst[i];
            il = nlast[i];
            double sum = 0.0;
            for (j = if_; j <= il; ++j) sum += pb[(j * (j + 1)) / 2];
            k = nat[i];
            pop[i] = sum;
            chrg[i] = tore[k] - pop[i];
        }
        std::fprintf(stdout, "\n\n\n         MULLIKEN POPULATIONS AND CHARGES\n\n");
        std::fprintf(stdout, "     NO.  ATOM   POPULATION      CHARGE\n");
        for (j = 1; j <= numat; ++j)
            std::fprintf(stdout, "%5d%3s%16.6f%14.6f\n", j,
                         elemnt[nat[j]].c_str(), pop[j], chrg[j]);
        if (icok > -1) {
            for (j = 1; j <= numat; ++j) std::fprintf(isyb, "%12.6f", pop[j]);
            for (j = 1; j <= numat; ++j) std::fprintf(isyb, "%12.6f", chrg[j]);
            std::fprintf(isyb, "\n");
        }
        to_screen("To_file: Mulliken");
    }
}

void mpcsyb(double* chr, int kchrge, double eionis, double& dip) {
    // F90: open(unit=isyb, file=syb_fn); writes actually go to unit 16.
    FILE* isyb = std::fopen(syb_fn.c_str(), "w");
    if (!isyb) { std::fprintf(stdout, "Error writing SYBYL MOPAC output\n"); return; }
    std::fprintf(isyb, "%4d%4d\n", 1, numat);
    for (int i = 1; i <= numat; ++i)
        std::fprintf(isyb, "%12.6f%12.6f%12.6f%12.6f\n",
                     coord[0][i], coord[1][i], coord[2][i], chr[i]);
    int nfilled = (nclose > nalpha) ? nclose : (nalpha > nbeta ? nalpha : nbeta);
    int i1 = (nfilled - 1 > 1) ? nfilled - 1 : 1;
    int i2 = (nfilled + 2 < norbs) ? nfilled + 2 : norbs;
    // F90 writes eigs(i1..i1+3): if i2 limits the upper window, keep the
    // F90 window starting at i1 (i1..i1+3); guarded by array size.
    int ehi = (i1 + 3 <= (int)eigs.size() - 1) ? i1 + 3 : (int)eigs.size() - 1;
    std::fprintf(isyb, "%12.6f%12.6f%12.6f%12.6f  %4d  HOMOs,LUMOs,# of occupied MOs\n",
                 eigs[i1], eigs[i1+1], eigs[i1+2], eigs[i1+3], nfilled);
    std::fprintf(isyb, "%12.6f%12.6f    HF and IP\n", escf, eionis);
    if (kchrge != 0) dip = 0.0;
    std::fprintf(isyb, "%4d%10.3f  Charge,Dipole Moment\n", kchrge, dip);
    if (keywrd.find(" MULL") != std::string::npos) mpcpop(1, isyb);
    else mpcpop(0, isyb);
    std::fclose(isyb);
}
