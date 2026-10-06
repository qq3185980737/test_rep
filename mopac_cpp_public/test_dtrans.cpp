// test_dtrans.cpp
#include <cstdio>
#include <vector>
#include "dtrans.h"
#include "symmetry_C.h"
namespace symmetry_C {
    int nclass=1;
    double elem[4][4][21] = {};
}
namespace funcon_C { double pi = 3.14159265358979; }
int main() {
    std::vector<double> d(6,0.0); d[3]=1.0;
    std::vector<std::vector<double>> r(3,std::vector<double>(3,0.0));
    for(int i=0;i<3;++i) r[i][i]=1.0;
    bool first=true;
    dtrans(d,1,first,r);
    bool ok=!first;
    std::printf("dtrans(skeleton) %s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}
