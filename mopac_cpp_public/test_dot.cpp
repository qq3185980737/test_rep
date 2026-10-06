// test_dot.cpp
#include <cstdio>
#include <vector>
#include "dot.h"
int main() {
    // Fortran 1-based: padding at index 0, data at 1..n.
    std::vector<double> x={0,1,2,3}, y={0,4,5,6};
    double s=dot(x,y,3);
    bool ok=(s==32.0);
    std::printf("dot=%g %s\n",s,ok?"PASS":"FAIL");
    return ok?0:1;
}
