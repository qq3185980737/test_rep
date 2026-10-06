// test_dtran2.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dtran2.h"
namespace funcon_C { double pi = 3.14159265358979; }
int main() {
    std::vector<std::vector<double>> r(3, std::vector<double>(3,0.0));
    for(int i=0;i<3;++i) r[i][i]=1.0;
    std::vector<std::vector<std::vector<double>>> t(13,
        std::vector<std::vector<double>>(5, std::vector<double>(5,0.0)));
    dtran2(r,t,1);
    double diag = t[1][2][2];
    bool ok = std::fabs(diag-1.0)<1e-9;
    std::printf("dtran2 identity T33=%g %s\n",diag, ok?"PASS":"FAIL");
    return ok?0:1;
}
