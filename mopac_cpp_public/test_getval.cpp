// test_getval.cpp
#include <cstdio>
#include <string>
#include "getval.h"
int main() {
    double x; std::string t;
    getval("1.234 H", x, t); bool a=(x==1.234 && t==" ");
    getval("R1 2.0", x, t); bool b=(x==-999.0 && t=="R1");
    std::printf("num x=%g sym t=%s %s\n",1.234,t.c_str(),(a&&b)?"PASS":"FAIL");
    return (a&&b)?0:1;
}
