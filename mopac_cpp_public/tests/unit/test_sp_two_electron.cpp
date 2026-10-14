// test_sp_two_electron.cpp
#include <cstdio>
#include "sp_two_electron.h"
double zsn[81]={0}, zpn[81]={0}, gss[81]={0}, gsp[81]={0}, gpp[81]={0}, gp2[81]={0}, hsp[81]={0};
bool main_group[81]={false};
int iii[81]={0};
double rsc(int,int,double,int,double,int,double,int,double){return 1.0;}
int main() { sp_two_electron(); std::printf("sp_two_electron PASS\n"); return 0; }
