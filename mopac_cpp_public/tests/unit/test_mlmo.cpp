// test_mlmo.cpp
#include "mlmo.h"
#include "MOZYME_C.h"
#include "molkst_C.h"
#include <cassert>
#include <cstdio>
#include <vector>
int main() {
    using namespace MOZYME_C;
    using namespace molkst_C;
    numat = 2; norbs = 4;
    iorbs.assign(norbs + 1, 0);
    iorbs[1] = 1; iorbs[2] = 4;
    ipad2 = 100; ipad4 = 100;
    std::vector<int> iz(numat + 1, 2), ib(numat + 1, 2);
    std::vector<int> nce(10, 0), ncf(10, 0);
    std::vector<int> ncocc(10, 0), ncvir(10, 0);
    std::vector<int> icocc(100, 0), icvir(100, 0);
    std::vector<double> cocc(200, 1.0), cvir(200, 1.0);
    int locc = 0, lvir = 0, nf_loc = 0, ne = 0, nocc = 0, nvir = 0;
    mlmo(locc, lvir, 1, 2, nf_loc, ne, nocc, nvir, iz.data(), ib.data(),
         nce.data(), ncf.data(), ncocc.data(), ncvir.data(), iorbs.data(),
         icocc.data(), icvir.data(), cocc.data(), cvir.data());
    assert(nocc == 1); assert(nvir == 1);
    assert(ncocc[1] == 0); assert(ncvir[1] == 0);
    assert(locc == 8);
    assert(lvir == 8);
    printf("mlmo diatomic PASS\n");
    return 0;
}