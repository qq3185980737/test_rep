// test_fock1.cpp
#include <cstdio>
#include "fock1.h"
#include "common_arrays_C.h"
#include "molkst_C.h"
int main() {
    using namespace common_arrays_C; using namespace molkst_C;
    numat=1;
    nfirst.assign(2,1); nlast.assign(2,1); nat.assign(2,6);
    std::vector<double> Fm={0},pt={0},Paa={0},Pbb={0};
    fock1(Fm,pt,Paa,Pbb);
    std::printf("fock1(s-atom) PASS\n");
    return 0;
}
