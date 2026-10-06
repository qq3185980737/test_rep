// newflg.cpp
#include "newflg.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
using common_arrays_C::na;
using common_arrays_C::nb;
using common_arrays_C::nc;
using common_arrays_C::coord;
using common_arrays_C::geo;
using common_arrays_C::txtatm;
using common_arrays_C::loc;
using common_arrays_C::xparam;
using molkst_C::natoms;
using molkst_C::numat;
using molkst_C::nvar;
extern void mopend(const char*);
#include <cstdio>

void newflg() {
    if (numat != natoms) {
        std::printf(" NEWGEO CAN ONLY BE USED IF THERE ARE NO DUMMY ATOMS\n");
        mopend("NEWGEO cannot be used here");
        return;
    }
    int found = 0;
    for (int i = 1; i <= numat; ++i) if (na[i] != 0) { found = 1; break; }
    if (!found) {
        std::printf(" There is a bug in MOZYME: NEWGEO cannot be used here\n");
        mopend("There is a bug in MOZYME");
        return;
    }
    int k = 0;
    for (int i = 1; i <= numat; ++i) {
        const std::string& s = txtatm[i];
        if (s.size() < 15) continue;
        if (s.substr(0,6) != "ATOM  ") continue;
        if (s.substr(12,3) != " N " && s.substr(12,3) != " C ") continue;
        na[i] = 0; nb[i] = 0; nc[i] = 0;
        for (int j = 1; j <= 3; ++j) geo[j][i] = coord[j-1][i];
        k = 1;
    }
    for (int i = 1; i <= nvar; ++i) xparam[i] = geo[loc[2][i]][loc[1][i]];
    if (k == 0) std::printf(" WARNING! NO BACKBONE ATOMS IDENTIFIED\n");
}
