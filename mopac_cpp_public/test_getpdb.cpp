// test_getpdb.cpp
#include <cstdio>
#include <vector>
#include "getpdb.h"
int main() { std::vector<std::vector<double>> g(4,std::vector<double>(1)); getpdb(g); std::printf("getpdb links OK PASS\n"); return 0; }
