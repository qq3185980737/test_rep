// test_setupk.cpp
#include <cstdio>
#include <vector>
#include "setupk.h"
namespace mozyme_C { extern std::vector<int> icocc,ncf,nncf,kopt; }
int main() {
    using namespace mozyme_C;
    kopt.assign(5,0); ncf.assign(3,0); nncf.assign(3,0);
    ncf[1]=2; nncf[1]=0; icocc.assign(3,0); icocc[1]=2; icocc[2]=3;
    setupk(1);
    bool ok=(kopt[2]==1 && kopt[3]==1 && kopt[1]==0);
    std::printf("kopt=(%d,%d,%d) %s\n",kopt[1],kopt[2],kopt[3],ok?"PASS":"FAIL");
    return ok?0:1;
}
