// test_linpack.cpp
#include <cstdio>
#include "linpack.h"
int main() { double a[16],det[2],w[4]; int ip[4],info; dgefa(a,4,2,ip,info); dgedi(a,4,2,ip,det,w,1); std::printf("linpack info=%d PASS\n",info); return 0; }
