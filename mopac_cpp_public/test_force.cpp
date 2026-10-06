// test_force.cpp
#include <cstdio>
#include <vector>
#include "force.h"
int main() {
    int n=2;
    std::vector<std::vector<double>> v(3,std::vector<double>(3,0.0));
    v[1][1]=-2; v[2][1]=1;
    v[1][2]=3; v[2][2]=-1;
    phase_lock(v,n);
    bool ok=(v[1][1]==2.0 && v[1][2]==3.0 && v[2][2]==-1.0);
    std::printf("phase_lock v11=%g v22=%g %s\n",v[1][1],v[2][2], ok?"PASS":"FAIL");
    return ok?0:1;
}
