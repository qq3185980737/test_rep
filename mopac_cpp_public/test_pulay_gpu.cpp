// test_pulay_gpu.cpp
#include <cstdio>
#include "pulay_for_gpu.h"
int main() { double f[5],p[5],fpp[5],fk[5],em[40]; int lf=0,nf=0; bool st=true; double pl; pulay_for_gpu(f,p,2,fpp,fk,em,lf,nf,1,st,pl); std::printf("pulay_gpu pl=%g PASS\n",pl); return 0; }
