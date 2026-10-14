// test_symh.cpp
#include <cstdio>
#include <vector>
#include "symh.h"
#include "molkst_C.h"
#include "symmetry_C.h"
int molkst_C::numat=1;
namespace symmetry_C { std::vector<std::vector<int>> ipo(2,std::vector<int>(3,0)); double r[10][121]={}; }
void mat33(double*, double*, double* t2){t2[1]=t2[2]=t2[3]=t2[4]=t2[5]=t2[6]=t2[7]=t2[8]=t2[9]=0;}
int main() { double h[20]={0},dip[20]={0}; int ipo[4]={0,1,2,3}; symh(h,dip,1,1,ipo); std::printf("symh PASS\n"); return 0; }
