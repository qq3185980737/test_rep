// test_symr.cpp
#include <cstdio>
#include <vector>
#include "symr.h"
#include "molkst_C.h"
int molkst_C::numat=1;
namespace common_arrays_C { std::vector<std::vector<double>> coord = {{0},{0},{0}}; }
int iw=6;
int main() { symr(); std::printf("symr PASS\n"); return 0; }
