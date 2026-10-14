// test_symtry.cpp
#include <cstdio>
#include "symtry.h"
int ndep=0;
int idepfn[4]={0}, locpar[4]={0}, locdep[4]={0}, na[4]={0};
double depmul[4]={0}, geo[12]={0};
void mopend(const char*){}
int main() { symtry(); std::printf("symtry PASS\n"); return 0; }
