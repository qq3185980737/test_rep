// test_superd.cpp
#include <cstdio>
#include "superd.h"
int natorb[20]={0,0,0,0,0,0,4};
const char* elemnt(int z){ return z==6?" C":" X"; }
int iw=6;
int main() { double c[4]={1,0,0,1},e[2]={-1,1}; int nat[1]={6}; superd(c,e,2,2,1,nat); std::printf("superd PASS\n"); return 0; }
