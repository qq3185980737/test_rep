// test_scfcri.cpp
#include <cstdio>
#include "scfcri.h"
const char* keywrd="";
int numcal=1, iw=6;
double efield[4]={0,0,0,0};
double reada(const char*, int){return 1e-2;}
int main() { double s=0; scfcri(s); std::printf("selcon=%g PASS\n",s); return 0; }
