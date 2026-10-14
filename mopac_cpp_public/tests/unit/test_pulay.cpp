// test_pulay.cpp
#include <cstdio>
#include "pulay.h"
#include "molkst_C.h"
void mamult(const double*,const double*,double*,int,double){}
void osinv(double*,int,double&d){d=1.0;}
using namespace molkst_C;
int main() {
    numcal=1; mpack=10;
    double f[10]={1,0,0,1,0,0,0,0,0,0},p[10]={1,0,0,1,0,0,0,0,0,0};
    double fpp[20],fk[20],em[400]={0};
    int lf=1,nf=1; bool st=true; double pl;
    pulay(f,p,2,fpp,fk,em,lf,nf,20,st,pl);
    std::printf("pulay pl=%g PASS\n",pl);
    return 0;
}
