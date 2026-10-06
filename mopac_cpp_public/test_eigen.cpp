// test_eigen.cpp
#include <cstdio>
#include "eigen.h"
#include "molkst_C.h"
int main() {
    using namespace molkst_C;
    norbs=10; nelecs=6; keywrd="";
    int no=0,nv=0;
    eigen_limits(no,nv);
    bool ok=(no==3 && nv==7);
    std::printf("eigen_limits no=%d nv=%d %s\n",no,nv, ok?"PASS":"FAIL");
    return ok?0:1;
}
