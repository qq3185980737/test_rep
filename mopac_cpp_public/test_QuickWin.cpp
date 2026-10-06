// test_QuickWin.cpp
#include <cstdio>
#include "QuickWin_bits.h"
int main() { double f[5],p[5],fpp[5],fk[5],em[40]; int lf=0,nf=0; bool st=true; double pl; pulay_for_gpu_qw(f,p,2,fpp,fk,em,lf,nf,1,st,pl); std::printf("qw pl=%g PASS\n",pl); return 0; }
