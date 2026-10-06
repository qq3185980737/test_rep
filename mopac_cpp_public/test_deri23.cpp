// test_deri23.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "deri23.h"
#include "meci_C.h"
#include "molkst_C.h"
namespace meci_C {
    int nmos, nelec;
    int nbo[4];
}
namespace molkst_C {
    int norbs, nopen;
    double fract;
}
int main() {
    using namespace meci_C;
    using namespace molkst_C;
    nbo[1]=1; nbo[2]=0; nbo[3]=1;
    nelec=0; nmos=1; norbs=2; fract=1.0;
    std::vector<double> e={0,1,3}, f={0,2.0}, fd(3,0), fci(3,0);
    std::vector<std::vector<double>> cmo(3, std::vector<double>(3,0));
    std::vector<double> emo(4,0);
    deri23(f, fd, e, fci, cmo, emo, 2, 1);
    bool ok = (std::abs(cmo[2][1]+1.0)<1e-9 && std::abs(cmo[1][2]-1.0)<1e-9);
    std::printf("cmo[2][1]=%g cmo[1][2]=%g %s\n", cmo[2][1], cmo[1][2], ok?"PASS":"FAIL");
    return ok?0:1;
}
