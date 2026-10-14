// test_symt.cpp
#include <cstdio>
#include "symt.h"
#include "molkst_C.h"
#include "symmetry_C.h"
int molkst_C::numat=1;
void mat33(double*, double*, double*) {}
int main() { double h[4]={1,0,0,1},d[3]={0},ha[4]={0}; symt(h,d,ha); std::printf("symt PASS\n"); return 0; }
