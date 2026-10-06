// test_diag.cpp
#include <cmath>
#include <cstdio>
#include <vector>
#include "diag.h"
int main() {
    std::vector<double> fao={0, 0.0,1.0,0.0};  // f11,f12,f22
    std::vector<std::vector<double>> v(3, std::vector<double>(3,0));
    v[1][1]=1; v[2][2]=1;
    std::vector<double> eig={0,-1.0,1.0};
    diag(fao, v, 1, eig, 3, 2);
    bool changed = (std::abs(v[1][1]-1.0)>1e-6);
    std::printf("v11=%g changed=%d %s\n", v[1][1], (int)changed, changed?"PASS":"FAIL");
    return changed?0:1;
}
