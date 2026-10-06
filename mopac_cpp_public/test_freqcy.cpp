// test_freqcy.cpp
#include <cstdio>
#include <vector>
#include "freqcy.h"
#include "to_screen_C.h"
namespace to_screen_C { std::vector<std::vector<double>> redmas; }
extern "C" void symt_(double*, double*, double*) {}
extern "C" void frame_(std::vector<double>&, int, int) {}
extern "C" void rsp_(double*, int, double*, double*) {}
extern "C" void phase_lock_(double*, int) {}
extern "C" void symtrz_(double*, double*, int, int) {}
int main() {
    int nvar=2;
    std::vector<double> f={0,1.0,2.0,3.0}, old={0,0,0,0};
    std::vector<double> w={0,0.5,2.0};
    freqcy_mass_weight(f,old,w,nvar);
    bool ok=(f[2]==2.0*0.5*2.0 && old[2]==2.0*1e5);
    std::printf("f2=%g old2=%g %s\n",f[2],old[2], ok?"PASS":"FAIL");
    return ok?0:1;
}
