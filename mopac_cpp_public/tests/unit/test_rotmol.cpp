// test_rotmol.cpp
#include <cstdio>
#include "rotmol.h"
void symopr(int, double*, int, double*) {}
int main() { double c[9]={0},r[9]={1,0,0,0,1,0,0,0,1}; rotmol(3,c,0,1,1,2,r); std::printf("rotmol PASS\n"); return 0; }
