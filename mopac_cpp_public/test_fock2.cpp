// test_fock2.cpp
#include <cstdio>
#include <vector>
#include "fock2.h"
int main() {
    int kr=0, ia=1, ib=1, ilim=1;
    std::vector<double> f={0,0}, pt={0,1}, pa={0,0.5};
    std::vector<std::vector<double>> w(2,std::vector<double>(2,0.0));
    w[1][1]=1.0;
    fock1dorbs(f,pt,pa,w,kr,ia,ib,ilim);
    bool ok=(f[1]==0.5 && kr==1);
    std::printf("fock1dorbs f=%g kr=%d %s\n",f[1],kr, ok?"PASS":"FAIL");
    return ok?0:1;
}
