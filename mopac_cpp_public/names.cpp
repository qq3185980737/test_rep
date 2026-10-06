// names.cpp — C++ translation of "names.F90".
// MOZYME residue/side-chain naming: walk backbone via nxtmer, tag atoms,
// compute phi/psi/omega dihedrals.

#include "names.h"

#include <cmath>
#include <string>
#include <vector>

#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "funcon_C.h"
#include "molkst_C.h"

// External helpers (ported elsewhere or stubbed).
void nxtmer(int iatom, int* nbackb);
void atomrs(int* lused, bool* ioptl, int ires, int n1, int io,
            int uni_res, bool first_res);
void dihed(const std::vector<std::vector<double>>& coord, int a, int b,
           int c, int d, double& angle);
void mopend(const std::string& s);

using namespace common_arrays_C;
using namespace MOZYME_C;
using namespace molkst_C;

void names(bool* ioptl, int* lused, int n1, int& ires, int nfrag,
           int io, int& uni_res, int& mres) {
    bool lreseq = keywrd.find("RESEQ") != std::string::npos;
    for (int i = 1; i <= 4; ++i) {
        nbackb[i] = 0;
        for (int j = 1; j <= 1; ++j) nxeno[i][j] = 0;
    }
    mxeno = 1;

    iatom = n1;
    int iold = -999;
    if (nfrag == 1)
        for (int i = 1; i <= numat; ++i) txtatm[i] = " ";
    int irold = ires + 1;
    bbone[1][ires + 1] = n1;
    ++ires;
    int ires_start = ires;
    bool first_res = true;
    std::vector<bool> iatoms(numat + 1, false);

    for (int ires_loop = ires_start; ires_loop <= maxres; ++ires_loop) {
        if (iatom == 0) break;  // goto 999
        if (iatoms[iatom]) goto loop_1000;
        iatoms[iatom] = true;

        ioptl[iatom] = true;
        nxtmer(iatom, nbackb.data());
        jatom = nbackb[3];

        ++uni_res;
        int i = io;
        if (!first_res) {
            i = 1;
            while (i <= numat && nat[i] == 8) ++i;
        }
        atomrs(lused, ioptl, ires, n1, i, uni_res, first_res);
        first_res = false;
        bbone[1][ires + 1] = jatom;
        bbone[2][ires] = nbackb[0];
        bbone[3][ires] = nbackb[1];

        if (loop != 1 && k > 0 && k != 23 && i != 1) {
            size_t sp = txeno[loop - 1].find(' ');
            int ii = (sp == std::string::npos ? 1 : (int)sp) - 1;
            if (ii != 1 && keywrd.find("ADD-H") == std::string::npos) {
                // "Residue: ires (XXX) contains a xeno group."
            }
            if (allres[ires].size() >= 4 && allres[ires][3] == '+') {
                // positive charge message
            } else if (allres[ires].size() >= 4 && allres[ires][3] == '-') {
                // negative charge message
            }
            allres[ires] = allres[ires].substr(0, 3) + "*";
        }
        if (iatom == jatom) goto loop_1000;
        if (iold == jatom) {
            mopend("Structure Unrecognizable");
            break;
        }
        iold = iatom;
        iatom = jatom;
        if (iatom == 0 && ires_loop == ires_start) goto loop_1010;
        ++ires;
    }
    return;

loop_1000:
    // Stray N at end.
    if (jatom > 0) {
        if (nat[jatom] == 7 && txtatm[jatom] == " ") {
            int l = 0;
            for (int i = 1; i <= nbonds[jatom]; ++i) {
                int kk = ibonds[i][jatom];
                if (!ioptl[kk] && nat[kk] > 1) ++l;
            }
            if (l == 0) {
                txtatm[jatom] = "ATOM  " + std::to_string(jatom) + "  N   UNK " +
                                std::to_string(ires + 1);
            }
        }
    }
    if (keywrd.find("ADD-H") != std::string::npos) return;

    if (ires > 0) {
        bool fr = true;
        for (int i = irold; i <= ires; ++i) {
            ++mres;
            if (!fr) {
                if (bbone[3][i - 1] * bbone[1][i] != 0 &&
                    bbone[2][i] * bbone[3][i] != 0)
                    dihed(coord, bbone[3][i - 1], bbone[1][i], bbone[2][i],
                          bbone[3][i], angles[1][mres]);
                else
                    angles[1][mres] = 0.0;
            }
            if (i != ires) {
                if (bbone[1][i] * bbone[2][i] != 0 &&
                    bbone[3][i] * bbone[1][i + 1] != 0)
                    dihed(coord, bbone[1][i], bbone[2][i], bbone[3][i],
                          bbone[1][i + 1], angles[2][mres]);
                else
                    angles[2][mres] = 0.0;
                if (bbone[2][i] * bbone[3][i] != 0 &&
                    bbone[2][i + 1] * bbone[1][i + 1] != 0)
                    dihed(coord, bbone[2][i], bbone[3][i], bbone[1][i + 1],
                          bbone[2][i + 1], angles[3][mres]);
                else
                    angles[3][mres] = 0.0;
            }
            for (int j = 1; j <= 3; ++j) {
                if (angles[j][mres] > funcon_C::pi)
                    angles[j][mres] -= 2.0 * funcon_C::pi;
                if (angles[j][mres] < -funcon_C::pi)
                    angles[j][mres] = 0.0;
                angles[j][mres] *= 180.0 / funcon_C::pi;
            }
            fr = false;
        }
        for (int i = 1; i <= ires; ++i) {
            int j;
            for (j = 1; j <= 20; ++j)
                if (allres[i].substr(0, 3) == afn[j].substr(0, 3)) break;
            if (j > 20) allr[i] = "X";
        }
    }

    // Label hydrogens.
    int nbreaks = 1;
    for (int j = 1; j <= numat; ++j) {
        if (nat[j] == 1) {
            if (j != breaks[nbreaks]) {
                if (ibonds[1][j] > 0) lused[j] = lused[ibonds[1][j]];
            } else {
                ++nbreaks;
            }
        }
    }
    if (numat != natoms + id && !lreseq) {
        int j = natoms - numat;
        int l = numat;
        for (int i = natoms - id; i >= 1; --i) {
            if (labels[i] != 99) {
                txtatm[i] = txtatm[l];
                lused[i] = lused[l];
                --l;
            } else {
                txtatm[i] = std::to_string(j);
                --j;
            }
        }
    }
    nres = ires;
    return;

loop_1010:
    nres = 1;
    if (allres[ires].substr(0, 3) != "   ") {
        // "AMINO ACID (residue: ires) = XXX"
    }
    for (int i = 1; i <= natoms; ++i)
        if (lused[i] != -1) lused[i] = 1;
    for (int j = 1; j <= numat; ++j)
        if (nat[j] == 1) lused[j] = lused[ibonds[1][j]];
    if (numat == natoms) return;
    {
        int j = natoms - numat;
        int l = numat;
        for (int i = natoms; i >= 1; --i) {
            if (labels[i] != 99) {
                txtatm[i] = txtatm[l];
                lused[i] = lused[l];
                --l;
            } else {
                txtatm[i] = std::to_string(j);
                --j;
            }
        }
    }
}
