// test_sort.cpp
#include <cstdio>
#include <complex>
#include "sort.h"
int main() {
    float v[5]={0,3,1,2,0};
    std::complex<float> vec[16];
    for (int i=0;i<4;++i) { vec[i*4]=std::complex<float>(v[i],0); vec[i*4+1]=std::complex<float>(v[i],1); vec[i*4+2]=std::complex<float>(v[i],2); vec[i*4+3]=std::complex<float>(v[i],3); }
    sort(v,vec,4);
    bool ok=(v[1]==0 && v[2]==1 && v[3]==2 && v[4]==3);
    std::printf("v=(%g,%g,%g,%g) %s\n",v[1],v[2],v[3],v[4],ok?"PASS":"FAIL");
    return ok?0:1;
}
