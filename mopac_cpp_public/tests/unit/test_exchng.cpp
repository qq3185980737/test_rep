// test_exchng.cpp
#include <cstdio>
#include <vector>
#include "exchng.h"
int main() {
    double b=0,d=0,q=0;
    std::vector<double> x={0,1,2}, y={0,0,0};
    exchng(3,b,4,d,5,q,x,y,2);
    bool ok=(b==3 && d==4 && q==5 && y[1]==1 && y[2]==2);
    std::printf("b=%g d=%g q=%g y2=%g %s\n",b,d,q,y[2], ok?"PASS":"FAIL");
    return ok?0:1;
}
