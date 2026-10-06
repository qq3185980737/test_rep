// test_genun.cpp
#include <cstdio>
#include <cmath>
#include <vector>
#include "genun.h"
int main() {
    int n=100;
    std::vector<std::vector<double>> u(4,std::vector<double>(n+1,0.0));
    genun(u,n);
    bool ok=(n>50 && n<=100);
    double nm=u[1][1]*u[1][1]+u[2][1]*u[2][1]+u[3][1]*u[3][1];
    ok &= (std::fabs(nm-1.0)<0.01);
    std::printf("genun n=%d norm=%g %s\n",n,nm, ok?"PASS":"FAIL");
    return ok?0:1;
}
