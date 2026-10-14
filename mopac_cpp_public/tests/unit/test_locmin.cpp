// test_locmin.cpp
#include <cstdio>
#include "locmin.h"
int main() { double x[5],p[5],ssq,alf,efs[5]; int nc=0; locmin(1,x,5,p,ssq,alf,efs,nc); std::printf("locmin links ssq=%g PASS\n",ssq); return 0; }
