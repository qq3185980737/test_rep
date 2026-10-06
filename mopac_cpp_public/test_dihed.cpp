// test_dihed.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "dihed.h"
namespace funcon_C { double pi = 3.14159265358979; }
int main() {
    std::vector<std::vector<double>> xyz(3, std::vector<double>(5,0.0));
    xyz[0][3]=1; xyz[2][2]=1; xyz[0][4]=-1;
    double ang;
    dihed(xyz, 3,2,1,4, ang);
    double deg = ang*180.0/3.14159265;
    bool ok = std::fabs(deg-180.0)<1.0;
    std::printf("dihedral=%.1f deg %s\n", deg, ok?"PASS":"FAIL");
    return ok?0:1;
}
