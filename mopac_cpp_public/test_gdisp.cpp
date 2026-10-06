// test_gdisp.cpp
#include <cstdio>
#include <cmath>
#include <vector>
#include "gdisp.h"
int main() {
    int natoms=2;
    std::vector<double> rcov(95,1.0);
    std::vector<int> nat={0,6,6};
    std::vector<std::vector<double>> xyz(4,std::vector<double>(3,0.0));
    xyz[1][1]=0; xyz[1][2]=2.0;
    std::vector<double> cn(3,0.0);
    ncoord(natoms,rcov,nat,xyz,cn);
    bool ok=(std::abs(cn[1]-0.5)<1e-9);
    std::printf("ncoord cn1=%g %s\n",cn[1], ok?"PASS":"FAIL");
    return ok?0:1;
}
