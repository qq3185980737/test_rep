// test_dofs.cpp
#include <cstdio>
#include <vector>
#include "dofs.h"
int main() {
    std::vector<std::vector<double>> eref(2, std::vector<double>(3,0.0));
    eref[1][1]=0.5; eref[1][2]=1.0;
    std::vector<double> dd(11,0.0);
    dofs(eref,1,2,dd,10,0.0,2.0);
    double tot=0; for(int i=1;i<=10;++i) tot+=dd[i];
    bool ok=(tot>0.0);
    std::printf("dofs total=%g %s\n",tot,ok?"PASS":"FAIL");
    return ok?0:1;
}
