// test_dftd3_bits.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dftd3_bits.h"
double eabh(int n,int A,int B,int H,const std::vector<std::vector<double>>& xyz,double sc,double cab);
int main() {
    // linear A-H-B, A=(-1,0,0), B=(1,0,0), H=(0,0,0): rab=2, rhm=0, cos=-1
    std::vector<std::vector<double>> xyz(4, std::vector<double>(4,0));
    xyz[1][1]=-1; xyz[1][3]=1; xyz[1][2]=0;
    double e = eabh(3, 1, 3, 2, xyz, 1.0, 1.0);
    int iadr=0,jadr=0;
    d3limit(250, 150, iadr, jadr);
    bool ok = (e < 0.0 && iadr==3 && jadr==2);
    std::printf("e=%g iadr=%d jadr=%d %s\n", e, iadr, jadr, ok?"PASS":"FAIL");
    return ok?0:1;
}
