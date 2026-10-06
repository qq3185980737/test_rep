// test_fock1_for_MOZYME.cpp
#include <cstdio>
#include <vector>
#include "fock1_for_MOZYME.h"
int main() {
    int iab=1, ilim=1, kr=0;
    std::vector<double> f={0,0}, pt={0,1};
    std::vector<std::vector<double>> w(2,std::vector<double>(2,0.0));
    w[1][1]=1.0;
    fock1_for_MOZYME(f,pt,w,kr,iab,ilim);
    bool ok=(f[1]==0.5 && kr==1);
    std::printf("fock1m f=%g kr=%d %s\n",f[1],kr, ok?"PASS":"FAIL");
    return ok?0:1;
}
