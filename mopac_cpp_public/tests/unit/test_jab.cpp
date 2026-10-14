// test_jab.cpp
#include <cstdio>
#include <cstring>
#include "jab.h"
int main() {
    double pja[17]={},pjb[17]={},w[101]={},f[200]={};
    pja[1]=1; pjb[1]=1; w[1]=1;
    jab(1,1,pja,pjb,w,f);
    // i=1 -> ioff=1, f[1]+=sumb1=1; joff=1 f[1]+=suma1=1 => f[1]=2
    bool ok=(f[1]==2.0);
    std::printf("jab f[1]=%g %s\n",f[1],ok?"PASS":"FAIL");
    return ok?0:1;
}
