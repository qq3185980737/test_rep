// test_xxx.cpp
#include <cstdio>
#include <string>
#include "xxx.h"
int main() {
    std::string r;
    xxx('R',1,2,3,0,r);
    bool ok=(r=="R123");
    std::printf("\"%s\" %s\n",r.c_str(),ok?"PASS":"FAIL");
    return ok?0:1;
}
