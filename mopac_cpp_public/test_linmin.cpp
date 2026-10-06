// test_linmin.cpp
#include <cstdio>
#include "linmin.h"
int main() { double x[10],p[10],a,f; bool ok; int ic; linmin(x,a,p,10,f,ok,ic,0.0); std::printf("linmin links a=%g PASS\n",a); return 0; }
