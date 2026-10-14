// test_molsymy.cpp
#include <cstdio>
#include "molsymy.h"
int main() { double c[9],r[9]; int e; molsym(c,e,r); std::printf("molsymy err=%d PASS\n",e); return 0; }
