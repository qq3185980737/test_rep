// test_kab.cpp
#include <cstdio>
#include "kab.h"
int main() {
    double pk[17]={},w[101]={},f[200]={0};
    pk[1]=1; w[1]=1;
    kab(1,2,pk,w,f);
    bool ok=(f[2]==-1.0);
    std::printf("kab f[2]=%g %s\n",f[2],ok?"PASS":"FAIL");
    return ok?0:1;
}
