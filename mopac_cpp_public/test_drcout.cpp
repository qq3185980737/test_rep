// test_drcout.cpp
#include <cstdio>
#include <vector>
#include "drcout.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    numcal=1; numat=1; keywrd="";
    std::vector<std::vector<double>> x(3,std::vector<double>(3,0)), g=x, v=x;
    std::vector<double> e3={0,10,1,0}, k3={0,0,1,0}, t3={0,10,2,0}, xt={0,0,0,0};
    std::vector<double> ch(1,0);
    int jl=0;
    drcout(x,g,v,3,0,e3,k3,t3,xt,1,ch,0.5,0,jl);
    std::printf("drcout(skeleton) PASS\n");
    return 0;
}
