// test_mbonds.cpp
#include "mbonds.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>
namespace molkst_C { int numat=0; int norbs=0; int nelecs=0; }
int ijbo(int i, int j) { return (i*(i-1))/2 + j - 1; }
void mopend(const char*) {}
int main() {
    using namespace molkst_C;
    numat = 2; norbs = 5; nelecs = 2;
    printf("start\n");
    std::vector<int> iorbs(numat+1, 0); iorbs[1]=1; iorbs[2]=4;
    std::vector<int> nfirst(numat+1,0); nfirst[1]=1; nfirst[2]=2;
    double catom[25] = {0};
    catom[0] = 1.0;
    catom[5] = 1.0;
    std::vector<double> f(20, 0.0);
    f[1] = 1.0; f[2] = -1.0; f[3] = 0.5;
    double cocc[100]={}, cvir[100]={};
    bool ui[5]={}, uj[5]={};
    bool lok=false;
    printf("calling mbonds\n");
    mbonds(0,0,f.data(),catom,nfirst.data(),1,2,ui,uj,lok,iorbs.data(),
           cocc,cvir,100,100,numat,norbs,5,100);
    printf("lok=%d\n",(int)lok);
    assert(lok);
    printf("cocc[1]=%f cvir[1]=%f\n",cocc[1],cvir[1]);
    printf("mbonds diatomic PASS\n");
    return 0;
}