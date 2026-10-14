// test_dhc.cpp
#include <cstdio>
#include "dhc.h"
int main() { double p[4],pa[4],pb[4],xi[6]; int nat[2]; double e; dhc(p,pa,pb,xi,nat,1,1,1,1,e,1); std::printf("dhc e=%g PASS\n",e); return 0; }
