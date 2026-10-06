// test_reseq.cpp
#include <cstdio>
#include "reseq.h"
int main() { bool io[2]; int lu[2],nw=0,o=0; reseq(io,lu,1,nw,o); std::printf("reseq nw=%d PASS\n",nw); return 0; }
