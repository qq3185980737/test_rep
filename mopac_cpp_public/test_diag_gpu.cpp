// test_dgpu.cpp
#include <cstdio>
#include "diag_for_GPU.h"
int main() { double fao[4],v[4],e[2]; diag_for_GPU(fao,v,1,e,2,4); std::printf("diag_gpu links OK PASS\n"); return 0; }
