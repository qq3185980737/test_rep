// test_alpb_and_xfac_mndod.cpp
#include <cmath>
#include <cstdio>
#include "alpb_and_xfac_mndod.h"
#include "parameters_C.h"
using namespace parameters_C;
static bool chk(const char* n, double g, double w){ bool ok=std::fabs(g-w)<1e-9; std::printf("%-10s %g vs %g %s\n",n,g,w,ok?"PASS":"FAIL"); return ok; }
int main(){
    alpb_and_xfac_mndod();
    bool ok=true;
    ok&=chk("alpb[11][1]",alpb[11][1],1.05225212);
    ok&=chk("alpb[13][13]",alpb[13][13],1.38788);
    ok&=chk("xfac[16][12]",xfac[16][12],1.0);
    ok&=chk("alpb[3][3]",alpb[3][3],0.0);
    return ok?0:1;
}
