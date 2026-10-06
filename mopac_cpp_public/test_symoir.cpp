// test_symoir.cpp
#include <cstdio>
#include "symoir.h"
#include <string>
#include <vector>
std::vector<int> jndex;
std::vector<std::string> namo;
int nirred=1, nclass=1, numat=1, lab=1;
std::string name="C1";
std::vector<std::string> jx={"","A","B"};
std::vector<std::string> group={"1","1"};
const char* keywrd="";
int iw=6;
double charmo(double*, int*, int, int, double*, int, bool){return 1;}
double charvi(double*, int, int, double*, int){return 1;}
double charst(double*, int*, int, int, double*, int, bool){return 1;}
int main() { jndex.assign(3,0);namo.assign(3,"A"); double v[4]={0},e[2]={0,0},r[9]={0}; symoir(1,v,e,2,r,2); std::printf("symoir PASS\n"); return 0; }
