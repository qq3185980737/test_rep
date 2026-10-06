// test_chrge_mz.cpp
#include <cstdio>
#include "chrge_for_MOZYME.h"
int ijbo(int i,int j){ return i*j-1; }
int main() {
    double p[10]={0,1,2,3,4,5,6,7,8,9};
    double q[2];
    int iorbs[3]={0,2,3};
    chrge_for_MOZYME(p,q,2,iorbs);
    bool ok=(q[1]==p[1]+p[3] && q[2]==p[4]+p[6]+p[9]);
    std::printf("q1=%g q2=%g %s\n",q[1],q[2],ok?"PASS":"FAIL");
    return ok?0:1;
}
