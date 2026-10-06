// test_alpb_and_xfac_pm3.cpp
#include <cmath>
#include <cstdio>
#include "alpb_and_xfac_pm3.h"
#include "parameters_C.h"
using namespace parameters_C;
static bool chk(const char* n, double g, double w){ bool ok=std::fabs(g-w)<1e-9; std::printf("%-10s %g vs %g %s\n",n,g,w,ok?"PASS":"FAIL"); return ok; }
int main(){
    alpb_and_xfac_pm3();
    bool ok=true;
    ok&=chk("alpb[11][1]",alpb[11][1],1.800472);
    ok&=chk("xfac[37][37]",xfac[37][37],2.654922);
    ok&=chk("alpb[56][56]",alpb[56][56],0.999997);
    ok&=chk("xfac[55][55]",xfac[55][55],8.549732);
    ok&=chk("alpb[6][6]",alpb[6][6],0.0);
    return ok?0:1;
}
