// test_dgpu.cpp
#include <cstdio>
#include "density_for_GPU.h"
int main() { double c[4],pp[4]; density_for_GPU(c,1,1,0,2,4,2,0,pp,0); std::printf("density_gpu links OK PASS\n"); return 0; }
