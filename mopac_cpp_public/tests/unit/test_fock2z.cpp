// test_fock2z.cpp
#include <cstdio>
#include <vector>
#include "fock2z.h"
int main() {
    int iab=1, jba=1, kr=0;
    std::vector<double> fii={0,0}, fjj={0,0}, fij={0,0};
    std::vector<double> pii={0,1}, pjj={0,2}, pij={0,0.5};
    std::vector<double> wj={0,1.0}, wk={0,1.0};
    focd2z(iab,jba,fii,fjj,fij,pii,pjj,pij,wj,wk,true,kr);
    bool ok=(fii[1]==2.0 && kr==1);
    std::printf("focd2z fii=%g kr=%d %s\n",fii[1],kr, ok?"PASS":"FAIL");
    return ok?0:1;
}
