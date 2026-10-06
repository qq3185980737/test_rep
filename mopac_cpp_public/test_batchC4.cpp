// test_batchC4.cpp — batch C4: lewis.F90 (lewis / remove_bond / check_h /
// check_CVS / txt_to_atom_no).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "lewis.h"
#include "molkst_C.h"
#include "mod_atomradii.h"
#include "mopend.h"

using namespace molkst_C;
using namespace common_arrays_C;
using namespace chanel_C;

// ---------- external dependencies (test-owned) ----------
double distance(int i, int j) {
    double dx = coord[1][i] - coord[1][j];
    double dy = coord[2][i] - coord[2][j];
    double dz = coord[3][i] - coord[3][j];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// set_up_dentate stub: builds raw connectivity from dentate_bonds.
std::vector<std::pair<int, int>> dentate_bonds;
void set_up_dentate() {
    for (auto& b : dentate_bonds) {
        nbonds[b.first] += 1;
        ibonds[nbonds[b.first]][b.first] = b.second;
        nbonds[b.second] += 1;
        ibonds[nbonds[b.second]][b.second] = b.first;
    }
}

static int cal_ctr = 0;
static void setup() {
    numcal = ++cal_ctr;
    moperr = false;
    mozyme = false;
    natoms = 0;
    numat = 0;
    maxtxt = 0;
    prt_topo = false;
    keywrd = " ";
    line = "";
    labels.assign(64, 0);
    nat.assign(64, 0);
    txtatm.assign(64, std::string(26, ' '));
    nbonds.assign(64, 0);
    ibonds.assign(16, std::vector<int>(64, 0));
    coord.assign(4, std::vector<double>(64, 0.0));
    mod_atomradii::radius.assign(64, 0.0);
    dentate_bonds.clear();
}

int main() {
    bool ok = true;

    // ---- T1: remove_bond drops the longest bond on both sides.
    {
        setup();
        numat = 4;
        nat[1] = 6; nat[2] = 6; nat[3] = 6; nat[4] = 6;
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 1; coord[2][2] = 0; coord[3][2] = 0;   // d=1
        coord[1][3] = 0; coord[2][3] = 2; coord[3][3] = 0;   // d=2
        coord[1][4] = 0; coord[2][4] = 0; coord[3][4] = 3;   // d=3 (longest)
        nbonds[1] = 3; ibonds[1][1] = 2; ibonds[2][1] = 3; ibonds[3][1] = 4;
        nbonds[2] = 1; ibonds[1][2] = 1;
        nbonds[3] = 1; ibonds[1][3] = 1;
        nbonds[4] = 1; ibonds[1][4] = 1;
        remove_bond(1);
        if (nbonds[1] != 2 || nbonds[4] != 0) {
            std::fprintf(stderr, "FAIL T1 nbonds %d %d\n", nbonds[1], nbonds[4]); ok = false;
        }
        if (nbonds[2] != 1 || ibonds[1][2] != 1) {
            std::fprintf(stderr, "FAIL T1 O2 %d %d\n", nbonds[2], ibonds[1][2]); ok = false;
        }
        if (ibonds[1][1] != 2 || ibonds[2][1] != 3) {
            std::fprintf(stderr, "FAIL T1 ibonds %d %d\n", ibonds[1][1], ibonds[2][1]); ok = false;
        }
    }

    // ---- T2: lewis nitro check bonds the two single-bonded O of N.
    {
        setup();
        numat = 3; natoms = 3;
        nat[1] = 7; nat[2] = 8; nat[3] = 8;   // N, O, O
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 1; coord[2][2] = 0; coord[3][2] = 0;
        coord[1][3] = -1; coord[2][3] = 0; coord[3][3] = 0;
        dentate_bonds = {{1, 2}, {1, 3}};   // N-O1, N-O2 only
        lewis(false);
        // N keeps 2; each O gains a second bond (O1-O2).  After the atomic-number
        // sort (descending), O(8) precedes N(7) in each bond list.
        if (nbonds[1] != 2 || nbonds[2] != 2 || nbonds[3] != 2) {
            std::fprintf(stderr, "FAIL T2 nbonds %d %d %d\n", nbonds[1], nbonds[2], nbonds[3]); ok = false;
        }
        if (ibonds[1][2] != 3 || ibonds[2][2] != 1 || ibonds[1][3] != 2 || ibonds[2][3] != 1) {
            std::fprintf(stderr, "FAIL T2 bridge %d %d %d %d\n",
                         ibonds[1][2], ibonds[2][2], ibonds[1][3], ibonds[2][3]); ok = false;
        }
    }

    // ---- T3a: check_h bonds an isolated H to nearest heavy atom, flags ibad.
    {
        setup();
        mozyme = true;
        numat = 3; natoms = 3;
        nat[1] = 1; nat[2] = 8; nat[3] = 6;
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 0; coord[2][2] = 0; coord[3][2] = 1.6;
        coord[1][3] = 0; coord[2][3] = 0; coord[3][3] = 3.0;
        mod_atomradii::radius[8] = 0.73;   // O covalent radius
        mod_atomradii::radius[6] = 0.77;
        int ibad = -1;
        check_h(ibad);
        if (nbonds[1] != 1 || ibonds[1][1] != 2) {
            std::fprintf(stderr, "FAIL T3a H bond %d %d\n", nbonds[1], ibonds[1][1]); ok = false;
        }
        if (nbonds[2] != 1 || ibonds[1][2] != 1) {
            std::fprintf(stderr, "FAIL T3a O %d %d\n", nbonds[2], ibonds[1][2]); ok = false;
        }
        if (ibad != 1) {
            std::fprintf(stderr, "FAIL T3a ibad=%d\n", ibad); ok = false;
        }
    }

    // ---- T3b: check_h removes the bridge bond of an H bonded to 2 O.
    {
        setup();
        mozyme = true;
        numat = 3; natoms = 3;
        nat[1] = 1; nat[2] = 8; nat[3] = 8;
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 0; coord[2][2] = 0; coord[3][2] = 1.2;
        coord[1][3] = 0; coord[2][3] = 0; coord[3][3] = 1.5;
        mod_atomradii::radius[8] = 0.73;
        nbonds[1] = 2; ibonds[1][1] = 2; ibonds[2][1] = 3;
        nbonds[2] = 1; ibonds[1][2] = 1;
        nbonds[3] = 1; ibonds[1][3] = 1;
        int ibad = -1;
        check_h(ibad);
        if (nbonds[1] != 1 || ibonds[1][1] != 2) {
            std::fprintf(stderr, "FAIL T3b H bond %d %d\n", nbonds[1], ibonds[1][1]); ok = false;
        }
        if (nbonds[3] != 0) {
            std::fprintf(stderr, "FAIL T3b O3 %d\n", nbonds[3]); ok = false;
        }
        if (ibad != 0) {
            std::fprintf(stderr, "FAIL T3b ibad=%d\n", ibad); ok = false;
        }
    }

    // ---- T4: check_CVS adds bond CVB(2,3).
    {
        setup();
        numat = 3; natoms = 3;
        nat[1] = 6; nat[2] = 6; nat[3] = 6;
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 1; coord[2][2] = 0; coord[3][2] = 0;
        coord[1][3] = 2; coord[2][3] = 0; coord[3][3] = 0;
        mod_atomradii::radius[6] = 0.77;
        keywrd = " CVB(2:3) ";
        check_CVS(false);
        if (nbonds[2] != 1 || ibonds[1][2] != 3 || nbonds[3] != 1 || ibonds[1][3] != 2) {
            std::fprintf(stderr, "FAIL T4 bond %d %d %d %d\n",
                         nbonds[2], ibonds[1][2], nbonds[3], ibonds[1][3]); ok = false;
        }
        if (moperr) {
            std::fprintf(stderr, "FAIL T4 moperr\n"); ok = false;
        }
    }

    // ---- T5: check_CVS deletes bond CVB(-2:-3); surviving atoms keep bonds,
    //          so no "all bonds deleted" fault.
    {
        setup();
        numat = 4; natoms = 4;
        nat[1] = 6; nat[2] = 6; nat[3] = 6; nat[4] = 6;
        coord[1][1] = 0; coord[2][1] = 0; coord[3][1] = 0;
        coord[1][2] = 1; coord[2][2] = 0; coord[3][2] = 0;
        coord[1][3] = 2; coord[2][3] = 0; coord[3][3] = 0;
        coord[1][4] = 3; coord[2][4] = 0; coord[3][4] = 0;
        mod_atomradii::radius[6] = 0.77;
        nbonds[2] = 2; ibonds[1][2] = 3; ibonds[2][2] = 4;
        nbonds[3] = 2; ibonds[1][3] = 2; ibonds[2][3] = 4;
        nbonds[4] = 2; ibonds[1][4] = 2; ibonds[2][4] = 3;
        keywrd = " CVB(-2:-3) ";
        check_CVS(false);
        if (nbonds[2] != 1 || ibonds[1][2] != 4 || nbonds[3] != 1 || ibonds[1][3] != 4) {
            std::fprintf(stderr, "FAIL T5 delete %d %d %d %d\n",
                         nbonds[2], ibonds[1][2], nbonds[3], ibonds[1][3]); ok = false;
        }
        if (moperr) {
            std::fprintf(stderr, "FAIL T5 moperr\n"); ok = false;
        }
    }

    // ---- T6: txt_to_atom_no resolves a PDB label to an atom number.
    {
        setup();
        maxtxt = 26;
        numat = 3; natoms = 3;
        txtatm[2] = std::string(12, ' ') + "CA1" + std::string(11, ' ');  // chars 13..15
        std::string text = "\"CA1\")";
        int m = -1;
        txt_to_atom_no(text, 1, false, m);
        if (text.find("2") == std::string::npos || text.find("CA1") != std::string::npos) {
            std::fprintf(stderr, "FAIL T6 text='%s'\n", text.c_str()); ok = false;
        }
    }

    if (ok) {
        std::fprintf(stderr, "ALL PASS\n");
        return 0;
    }
    std::fprintf(stderr, "FAILED\n");
    return 1;
}
