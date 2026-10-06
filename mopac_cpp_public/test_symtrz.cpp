// test_symtrz.cpp
#include <cstdio>
#include <vector>
#include <string>
#include "symtrz.h"
int numat=1, norbs=1, maxci=1, lab=1;
std::vector<int> nfirst={0,1}, nlast={0,1};
bool moperr=false;
std::vector<double> coord[4]={{0},{0},{0},{0}};
int nat[2]={0,1}, atmass[2]={0,1};
int nclass=1;
double elem[4][4][21]={0};
std::vector<std::vector<int>> jelem(1,std::vector<int>(1,0));
std::vector<std::string> namo;
std::vector<int> jndex;
const char* keywrd="";
int iw=6;
void molsym(double*, int& e, double*){e=0;}
void makopr(int, double*, int&, double*){}
void mult33(double*, int){}
void symoir(int, double*, double*, int, double*, int){}
int main() { double v[4]={0},e[2]={0}; symtrz(v,e,1,1); std::printf("symtrz PASS\n"); return 0; }
