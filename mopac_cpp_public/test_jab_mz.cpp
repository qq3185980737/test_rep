// test_jab_mz.cpp
#include <cstdio>
#include "jab_for_MOZYME.h"
int main() {
    double pja[17]={},pjb[17]={},w[101]={},f1[200]={},f2[200]={};
    pja[1]=1;pjb[1]=1;w[1]=1;
    jab_for_MOZYME(1,1,pja,pjb,w,f1,f2);
    bool ok=(f1[1]==1.0 && f2[1]==1.0);
    std::printf("f1=%g f2=%g %s\n",f1[1],f2[1],ok?"PASS":"FAIL");
    return ok?0:1;
}
