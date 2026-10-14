// test_getgeo.cpp
#include <cstdio>
#include <vector>
#include "getgeo.h"
int main() { std::vector<int> l(1),na(1),nb(1),nc(1); std::vector<std::vector<double>> g(4,std::vector<double>(1)),x(4,std::vector<double>(1)); std::vector<std::vector<int>> lo(4,std::vector<int>(1)); getgeo(5,l,g,x,lo,na,nb,nc,false); std::printf("getgeo links OK PASS\n"); return 0; }
