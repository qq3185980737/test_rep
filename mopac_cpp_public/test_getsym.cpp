// test_getsym.cpp
#include <cstdio>
#include <vector>
#include "getsym.h"
int main() { std::vector<int> a(100),b(100),c(100); std::vector<double> d(100); getsym(a,b,c,d); std::printf("getsym links OK PASS\n"); return 0; }
