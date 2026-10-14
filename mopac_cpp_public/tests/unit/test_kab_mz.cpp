// test_kab_mz.cpp
#include <cstdio>
#include "kab_for_MOZYME.h"
int main() {
    double pk[17]={},w[101]={},f[200]={0};
    pk[1]=1; w[1]=1;
    kab_for_MOZYME(1,1,pk,w,f);
    bool ok=(f[1]==-0.5);
    std::printf("kab_mz f[1]=%g %s\n",f[1],ok?"PASS":"FAIL");
    return ok?0:1;
}
