// test_mopac.cpp
#include <cstdio>
#include "mopac.h"
int main() { int r=mopac_main(0,nullptr); std::printf("mopac_main r=%d PASS\n",r); return 0; }
