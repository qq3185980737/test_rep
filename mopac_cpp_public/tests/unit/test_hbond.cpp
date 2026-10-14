// test_hbond.cpp — numeric checks for H_bond_correction_bits.
#include "H_bond_correction_bits.h"
#include "bangle.h"
#include "common_arrays_C.h"
#include "dihed.h"
#include "molkst_C.h"
#include "elemts_C.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using namespace common_arrays_C;
using namespace elemts_C;
using namespace molkst_C;

// minimal stub for reada (full version lives in reada.F90)
double reada(const std::string&, int) { return 0.0; }

int main() {
    int ok = 1;
    // truncation
    if (truncation(2.0, 4.0, 1.0) != 4.0) { std::printf("FAIL trunc a\n"); ok = 0; }
    if (std::fabs(truncation(4.0, 4.0, 1.0) - 4.25) > 1e-12) { std::printf("FAIL trunc b: %g\n", truncation(4.0,4.0,1.0)); ok = 0; }
    if (truncation(6.0, 4.0, 1.0) != 6.0) { std::printf("FAIL trunc c\n"); ok = 0; }
    // geometry setup: water dimer O1-H2...O3 along x
    id = 0;
    numat = 3;
    nat.assign(4, 0); nat[1]=8; nat[2]=1; nat[3]=8;
    coord.assign(4, std::vector<double>(4, 0.0));
    coord[1][1]=0.0; coord[2][1]=0.0; coord[3][1]=0.0;
    coord[1][2]=0.96; coord[2][2]=0.0; coord[3][2]=0.0;
    coord[1][3]=2.8; coord[2][3]=0.0; coord[3][3]=0.0;
    elemnt.assign(95, "  "); elemnt[1]=" H "; elemnt[8]=" O ";
    txtatm.assign(4, " ");
    // distance / connected
    double dd = distance_hb(1, 3);
    if (std::fabs(dd - 2.8) > 1e-12) { std::printf("FAIL distance: %g\n", dd); ok = 0; }
    if (!connected_hb(1, 2, 4.0) || std::fabs(Rab - 0.96) > 1e-12) { std::printf("FAIL connected\n"); ok = 0; }
    // find X-H bonds: O1-H2 within 1.15^2, O3-H2 = 1.84 > threshold -> only 1 H
    {
        std::vector<int> acc(8, 0), h_b(8, 0);
        int nacc = 0, nhb = 0;
        find_XH_bonds(acc, nacc, h_b, nhb);
        if (nacc != 2 || nhb != 1 || h_b[1] != 2) { std::printf("FAIL find_XH_bonds: nacc=%d nhb=%d\n", nacc, nhb); ok = 0; }
    }
    // all_h_bonds: exactly one H-bond triple (O1,H2,O3)
    {
        std::vector<int> l1(16,0), l2(16,0), l3(16,0);
        int nrpairs = 0;
        all_h_bonds(l1, l2, l3, 16, nrpairs);
        if (nrpairs != 1) { std::printf("FAIL all_h_bonds nrpairs=%d\n", nrpairs); ok = 0; }
        else if (!(l1[1]==1 && l2[1]==2 && l3[1]==3)) { std::printf("FAIL triple: %d %d %d\n", l1[1], l2[1], l3[1]); ok = 0; }
    }
    // prt_hbonds: energy < -0.5 registers a bond
    {
        maxtxt = 1;
        numcal = 5;
        keywrd = "";
        prt_hbonds(1, 2, 3, -1.0);
        if (P_Hbonds != 1 || H_energy[1] != -1.0 || H_txt[1].empty()) { std::printf("FAIL prt_hbonds\n"); ok = 0; }
    }
    // bonding
    {
        double covrad[95] = {0};
        covrad[1] = 0.32; covrad[8] = 0.71;
        if (std::fabs(bonding(1, 3, covrad) - 1.42) > 1e-12) { std::printf("FAIL bonding\n"); ok = 0; }
    }
    std::printf(ok ? "ALL PASS\n" : "FAILED\n");
    return ok ? 0 : 1;
}
