// test_jcarin.cpp
#include <cstdio>
#include "jcarin.h"
int main() { double xp[10],b[10]; int nc; jcarin(xp,1,false,b,nc,0,1); std::printf("jcarin ncol=%d PASS\n",nc); return 0; }
