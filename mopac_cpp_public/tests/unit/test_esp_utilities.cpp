// test_esp_utilities.cpp
#include <cstdio>
#include <vector>
#include "esp_utilities.h"
int main() {
    std::vector<float> sx={0,1,2,3}, sy={0,1,1,1};
    saxpy(3,2.0f,sx,1,sy,1);
    bool ok=(sy[1]==3 && sy[2]==5 && sy[3]==7);
    float d=sdot(3,sx,1,sx,1);
    ok &= (d==14.0f);
    float n=snrm2(3,sx,1);
    ok &= (n>3.7 && n<3.8);
    std::printf("saxpy=%g sdot=%g snrm2=%g %s\n",sy[3],d,n, ok?"PASS":"FAIL");
    return ok?0:1;
}
