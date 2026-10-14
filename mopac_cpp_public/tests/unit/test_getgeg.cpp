// test_getgeg.cpp
#include <cstdio>
#include <vector>
#include "getgeg.h"
int main() { std::vector<int> l(1),na(1),nb(1),nc(1); std::vector<std::vector<double>> g(4,std::vector<double>(1)); std::vector<std::vector<int>> lo(4,std::vector<int>(1)); getgeg(5,l,g,lo,na,nb,nc); std::printf("getgeg links OK PASS\n"); return 0; }
