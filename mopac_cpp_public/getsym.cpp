// getsym.cpp — C++ translation of MOPAC 2016 "getsym.F90".
// Reads parameter-dependence (symmetry) data and fills locpar/idepfn/locdep/
// depmul. Input lines are read from standard input (maps to Fortran unit ir).
#include "getsym.h"
#include "molkst_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

namespace {
using molkst_C::natoms;
using molkst_C::ndep;
using molkst_C::line;
using chanel_C::iw;
using common_arrays_C::na;
}  // namespace

void nuchar(char* line, int l_line, double* value, int& nvalue);

static const char* texti[19] = {
    " BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ",
    " BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ",
    " DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE",
    " DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS 270 DEGREES + REFERENCE DIHEDRAL  ",
    " DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ",
    " BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ",
    " BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ",
    " BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ",
    " DO NOT USE - USE SYMMETY FUNCTION 19 INSTEAD               ",
    " BOND LENGTH IS A MULTIPLE OF THE REFERENCE BOND LENGTH     "};
static const char* textx[19] = {
    " X COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ",
    " X COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ",
    " X COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ",
    " X COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ",
    " X COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ",
    " X COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ",
    " Y COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ",
    " Z COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ",
    " NOT USED                                                   "};

extern void mopend(const std::string&);

void getsym(std::vector<int>& locpar, std::vector<int>& idepfn,
            std::vector<int>& locdep, std::vector<double>& depmul) {
    int ivalue[101] = {0};
    double value[101] = {0.0};
    std::string used[39];
    int n_used = 0;
    int nerror = 0;
    int n = 0;

    std::fprintf(stdout, "\n\n\n%33s\n\n%12s%19s%16s\n",
                 "PARAMETER DEPENDENCE DATA", "REFERENCE ATOM", "FUNCTION NO.",
                 "DEPENDENT ATOM(S)");

    if (ndep > 0) {
        int j = 2;
        for (int i = 1; i <= ndep; ++i) {
            ++j;
            ivalue[j] = locdep[i];
            bool ok = (i == ndep);
            if (!ok) ok = (locpar[i] != locpar[i + 1] || idepfn[i] != idepfn[i + 1]);
            if (ok) {
                std::fprintf(stdout, "%13d%19d%16d%6d%6d%6d%6d%6d%6d", locpar[i], idepfn[i],
                             ivalue[3], ivalue[4], ivalue[5], ivalue[6], ivalue[7], ivalue[8]);
                for (int l = 9; l <= j; l += 7) {
                    std::fprintf(stdout, "\n%43s%7d%7d%7d%7d%7d%7d%7d", "",
                                 ivalue[l], ivalue[l + 1], ivalue[l + 2], ivalue[l + 3],
                                 ivalue[l + 4], ivalue[l + 5], ivalue[l + 6]);
                }
                std::fprintf(stdout, "\n");
                j = 2;
            }
        }
        goto done;
    }
    depmul[1] = 0.0;
    for (;;) {
        std::string buf;
        if (!std::getline(std::cin, buf)) goto done;
        line = buf;
        // Trim trailing blanks.
        size_t e = line.find_last_not_of(' ');
        if (e == std::string::npos) line = "";
        else line = line.substr(0, e + 1);
        int nvalue = 0;
        std::vector<char> cb(line.begin(), line.end());
        cb.push_back('\0');
        nuchar(cb.data(), (int)line.size(), value, nvalue);
        for (int i = 1; i <= nvalue; ++i) ivalue[i] = (int)std::lround(value[i]);
        if (nvalue == 0 || std::fabs(value[3]) < 1.e-20) goto done;
        if (ivalue[2] == 19) {
            if (na[ivalue[1]] == 0) {
                std::fprintf(stdout, "Atom %d is Cartesian.  Function 19 cannot be used here.\n", ivalue[1]);
                mopend("Error in Symmetry Data");
                return;
            }
            for (int i = 4; i <= nvalue; ++i) {
                if (ivalue[i] == 0) break;
                ++ndep;
                locdep[ndep] = ivalue[i];
                locpar[ndep] = ivalue[1];
                idepfn[ndep] = 19;
                ++n;
                double sum = value[3] * value[3];
                for (int l = 1; l <= 24; ++l) {
                    int j = (int)std::lround(sum * l);
                    if (std::fabs(j - sum * l) < l * 1.e-4) {
                        value[3] = std::sqrt((1.0 * j) / l);
                        break;
                    }
                }
                depmul[n] = value[3];
            }
        } else {
            if (na[ivalue[1]] != 0 && ivalue[2] == 18) {
                std::fprintf(stdout, "Atom %d is internal.  Function 18 cannot be used here.\n", ivalue[1]);
                mopend("Error in Symmetry Data");
                return;
            }
            for (int i = 3; i <= nvalue; ++i) {
                if (ivalue[i] == 0) break;
                ++ndep;
                locdep[ndep] = ivalue[i];
                locpar[ndep] = ivalue[1];
                idepfn[ndep] = ivalue[2];
                if (ivalue[i] > natoms) nerror = 1;
            }
        }
        int ll = nvalue;  // Fortran: i - 1 after loop (i = nvalue+1)
        while (ll > 0 && ivalue[ll] == 0) --ll;
        if (ivalue[2] == 19) {
            std::fprintf(stdout, "%13d%13d%13.8f%9d", ivalue[1], ivalue[2], value[3], ivalue[4]);
            for (int j = 5; j <= ll; ++j) std::fprintf(stdout, "%6d", ivalue[j]);
            std::fprintf(stdout, "\n");
        } else {
            std::fprintf(stdout, "%13d%19d%16d", ivalue[1], ivalue[2], ivalue[3]);
            for (int j = 4; j <= ll; ++j) std::fprintf(stdout, "%6d", ivalue[j]);
            std::fprintf(stdout, "\n");
        }
        int ii = (na[ivalue[1]] == 0) ? 2 : 1;
        std::string desc = (ii == 1) ? texti[ivalue[2] - 1] : textx[ivalue[2] - 1];
        int jj = 1;
        for (; jj <= n_used; ++jj)
            if (used[jj - 1] == desc) break;
        if (jj > n_used) {
            ++n_used;
            used[n_used - 1] = desc;
        }
        if (nerror == 1) {
            std::fprintf(stdout, " A SYMMETRY FUNCTION IS USED TO DEFINE A NON-EXISTENT ATOM\n");
            mopend("A SYMMETRY FUNCTION IS USED TO DEFINE A NON-EXISTENT ATOM");
            return;
        }
    }
done:
    std::fprintf(stdout, "\n%23s\n", "DESCRIPTIONS OF THE FUNCTIONS USED");
    for (int j = 1; j <= 18; ++j) {
        for (int i = 1; i <= n_used; ++i)
            if (used[i - 1] == texti[j - 1]) {
                std::fprintf(stdout, "%4d%5s%s\n", j, "", texti[j - 1]);
                break;
            }
    }
    for (int j = 1; j <= 18; ++j) {
        for (int i = 1; i <= n_used; ++i)
            if (used[i - 1] == textx[j - 1]) {
                std::fprintf(stdout, "%4d%5s%s\n", j, "", textx[j - 1]);
                break;
            }
    }
}

// No-arg stub: readmo.cpp calls getsym() without symmetry variables.
// TODO: real translation of symmetry.F90 constraint handling (M05 batch).
void getsym() {}
