// test_cublas.cpp
#include <cstdio>
#include "mod_calls_cublas.h"
int main() { double v[1]={1},r; call_asum_cublas(1,v,1,&r); std::printf("asum=%g PASS\n",r); return 0; }
