// test_set_up_rapid.cpp
#include <cstdio>
#include "set_up_rapid.h"
int numred=0, mode=0, nelred=0, numat=1, nelecs=2;
int step_num=0, iflepo=0;
bool use_ref_geo=false, moperr=false;
double escf=0, grad[3]={0}, xparam[3]={0};
const char* keywrd="";
void picopt(int){}
void compfg(double*, bool, double&, bool, double*, bool){}
void pinout(int){}
void hcore_for_MOZYME(){}
void dcart(double*, double*){}
int main() { set_up_rapid("OF"); set_up_rapid("RE"); std::printf("rapid PASS\n"); return 0; }
