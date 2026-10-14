// test_names.cpp
#include <cstdio>
#include "names.h"
int main() { bool io[2]; int lu[2],ires=0,ur=0,mr=0; names(io,lu,1,ires,1,1,ur,mr); std::printf("names ires=%d PASS\n",ires); return 0; }
