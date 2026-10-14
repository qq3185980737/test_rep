// test_swap.cpp
#include <cstdio>
#include "swap.h"
int iw=6;
int main() { double c[4]={1,0,0,1}; int ifill=-1; swap(c,2,2,1,ifill); std::printf("swap ifill=%d PASS\n",ifill); return 0; }
