// test_lbfgs.cpp
#include <cstdio>
#include "lbfgs.h"
int main() { double x[10],e; lbfgs(x,e); std::printf("lbfgs links e=%g PASS\n",e); return 0; }
