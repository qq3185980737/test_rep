// test_tidy.cpp
#include <cstdio>
#include "tidy.h"
int iorbs[10]={0,1,1}, jopt[10]={0,1,1};
double thresh=1e-8;
int numat=1, step_num=2, norbs=2, numred=1, numcal=1;
bool moperr=false;
const char* keywrd="";
int iw=6;
void selmos(int, int*, int*, int, double*, int, int*, int*, int, int, int*, int*, int*, int*, int*, int){}
void memory_error(const char*){}
void pinout(int){}
int main() { int nc[3]={0,1,1},ic[5]={0,1,2,0,0},n01=5,ln=0,mn=0,nnc[3]={0},ncmo[3]={0}; double c[5]={1,0,0,0,0}; tidy(2,nc,ic,n01,c,5,nnc,ncmo,ln,mn,1); std::printf("tidy ln=%d mn=%d PASS\n",ln,mn); return 0; }
