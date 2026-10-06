// test_compct.cpp
#include <cstdio>
#include <vector>

#include "compct.h"

int main() {
    int nmos = 3, idone = 3, itop = 3;
    std::vector<int> nncnew = {0,0,0,20}, ncmnew = {0,0,0,200};
    std::vector<int> ncnew(4, 0);
    std::vector<int> nc = {0, 2, 0, 3};
    std::vector<int> ic(30, 0);
    std::vector<int> iws = {0, 5, 0, 7};
    std::vector<double> c(300, 0.0);
    ic[21] = 11; ic[22] = 22;
    int LB = 0, MB = 0;
    compct(nncnew, ncnew, ncmnew, itop, nc, ic, iws, 30, c, 300, nmos, idone,
           LB, MB, 10, 100);
    bool ok = (LB == 8 && MB == 95 && nc[1] == 0 && nc[2] == 2 &&
               ic[9] == 11 && ic[10] == 22 && iws[2] == 5);
    std::printf("lb=%d mb=%d nc=(%d,%d) ic[9..10]=(%d,%d) %s\n",
                LB, MB, nc[1], nc[2], ic[9], ic[10], ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
