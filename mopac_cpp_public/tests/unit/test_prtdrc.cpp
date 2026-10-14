// test_prtdrc.cpp
#include <cstdio>
#include "prtdrc.h"
int main() { double xp[3],rf[3],v[3],gt=0,et=0; int mc[6]; prtdrc(0,xp,rf,0,gt,et,v,mc,1,0); std::printf("prtdrc links OK PASS\n"); return 0; }
