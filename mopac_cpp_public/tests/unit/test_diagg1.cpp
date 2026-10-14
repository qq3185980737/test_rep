// test_diagg1.cpp
#include <cstdio>
#include <vector>
#include "diagg1.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal=1; keywrd=""; numat=1;
    std::vector<double> fao(3), eigv(2), ws(2), fmo(2), avir(2), aocc(2), aov(2);
    std::vector<char> latoms(2);
    std::vector<std::vector<int>> ifmo(3, std::vector<int>(3,0));
    int nij=0;
    diagg1(fao, 0, 0, eigv, ws, latoms, ifmo, fmo, 2, nij, 1, avir, aocc, aov);
    bool ok = (aocc[0]==0.0);
    std::printf("diagg1(skeleton) %s\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
