// test_atomrs.cpp
#include <cmath>
#include <cstdio>
#include <vector>

#include "atomrs.h"
#include "common_arrays_C.h"

using namespace common_arrays_C;

int main() {
    // Build a graph:
    //   atom 3 = peptide N (nat=7), bonded to H(1), C-alpha(4), C=O(5)
    //   C-alpha(4) bonded to N(3), H(7), sidechain C(8)
    //   C=O carbon(5) bonded to N(3), O(9), prev-N(10)
    //   O(9) terminal (nbonds=1)
    const int M = 12;
    nat.assign(M + 1, 0);
    nbonds.assign(M + 1, 0);
    ibonds.assign(5, std::vector<int>(M + 1, 0));
    auto bond = [&](int a, int b) {
        ++nbonds[a];
        ibonds[nbonds[a]][a] = b;
    };
    nat[1] = 1; nat[3] = 7; nat[4] = 6; nat[5] = 6;
    nat[7] = 1; nat[8] = 6; nat[9] = 8; nat[10] = 7;

    bond(3, 1); bond(3, 4); bond(3, 5);     // N: H, Cα, C'
    bond(4, 3); bond(4, 7); bond(4, 8);     // Cα: N, H, Cβ
    bond(5, 3); bond(5, 9); bond(5, 10);    // C': N, O, prev-N
    bond(9, 5);                              // terminal O

    bool p3 = peptide_n(3);
    bool p4 = peptide_n(4);  // C, not N
    bool p1 = peptide_n(1);  // H
    std::printf("peptide_n(3)=%d (expect 1) %s\n", (int)p3, p3 ? "PASS" : "FAIL");
    std::printf("peptide_n(4)=%d (expect 0) %s\n", (int)p4, !p4 ? "PASS" : "FAIL");
    std::printf("peptide_n(1)=%d (expect 0) %s\n", (int)p1, !p1 ? "PASS" : "FAIL");

    // Negative case: N bonded to O-bearing C that also has another bond
    // (N-C-O-R) → should reject. Make O(9) have a second bond.
    bond(9, 8);
    bool p3b = peptide_n(3);
    std::printf("peptide_n(3) with O-R=%d (expect 0) %s\n", (int)p3b, !p3b ? "PASS" : "FAIL");

    return (p3 && !p4 && !p1 && !p3b) ? 0 : 1;
}
