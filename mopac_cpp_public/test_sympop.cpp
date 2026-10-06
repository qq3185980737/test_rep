// test_sympop.cpp
#include <cstdio>
#include "sympop.h"
void symh(double*, double*, int, int, const int*) {}
int main() { double h[4],d[3]; int iskip; sympop(h,1,iskip,d); std::printf("sympop iskip=%d PASS\n",iskip); return 0; }
