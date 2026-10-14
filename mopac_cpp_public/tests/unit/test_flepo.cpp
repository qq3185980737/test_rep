// test_flepo.cpp
#include <cstdio>
#include <vector>
#include "flepo.h"
static double ddot_(const std::vector<double>& a,const std::vector<double>& b){double s=0;for(int i=1;i<=2;i++)s+=a[i]*b[i];return s;}
int main() {
    int n=2;
    std::vector<double> xvar={0,1,1}, gvar={0,2,2}, gg={0,2,2};
    std::vector<double> H={0,0,0};
    double sy=ddot_(xvar,gvar);
    double yhy=ddot_(gg,gvar);
    flepo_update_H(n,xvar,gvar,gg,sy,yhy,false,H);
    bool ok=(H[3]==-0.25);
    std::printf("BFGS H22=%g sy=%g %s\n",H[3],sy, ok?"PASS":"FAIL");
    return ok?0:1;
}
