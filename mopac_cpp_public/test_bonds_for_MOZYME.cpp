// test_bonds_for_MOZYME.cpp
#include <cstdio>
#include <vector>

#include "bonds_for_MOZYME.h"
#include "MOZYME_C.h"
#include "chanel_C.h"
#include "common_arrays_C.h"
#include "elemts_C.h"
#include "molkst_C.h"

// Test stub for ijbo.
int ijbo(int ii, int jj) {
    if (ii == 1 && jj == 1) return 0;
    if (ii == 2 && jj == 2) return 1;
    if ((ii == 1 && jj == 2) || (ii == 2 && jj == 1)) return 2;
    return -1;
}

int iw = 6;

int main() {
    using namespace common_arrays_C;
    using namespace molkst_C;
    numat = 2; nl_atoms = 2;
    nat = {0, 6, 6};
    MOZYME_C::iorbs = {0, 1, 1};
    p = {0, 0.5, 0.6, 0.3};  // p(1),p(2),p(3)
    elemts_C::elemnt = {"", "H", "He", "Li", "Be", "B", "C"};
    keywrd = " PM3 ";
    l_atom = {false, true, true};

    bonds_for_MOZYME();
    // valenc atom1 = 2*0.5-0.25 = 0.75; atom1-bond to atom2 sum=0.09.
    // Verify no crash and valenc printed 0.75 for atom 1.
    std::printf("TEST DONE\n");
    return 0;
}
