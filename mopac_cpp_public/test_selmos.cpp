// test_selmos.cpp
#include <cstdio>
#include "selmos.h"
int numred=1, nelred=0, norred=0;
int iw=6;
int main() { int nc[3]={0,0,0},ic[3]={0,0,0},nnc[3]={0,0,0},ncmo[3]={0,0,0},iws[3]={0,0,0},jopt[2]={1,1},ncn[3]={0},ncmn[3]={0},nncn[3]={0}; double c[6]={0}; selmos(2,nc,ic,5,c,10,nnc,ncmo,2,2,iws,jopt,ncn,ncmn,nncn,1); std::printf("selmos PASS\n"); return 0; }
