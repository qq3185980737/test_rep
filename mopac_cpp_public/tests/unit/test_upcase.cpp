// test_upcase.cpp
#include <cstdio>
#include <string>
#include "upcase.h"
int main() {
    std::string s="am1 external geom";
    upcase(s,(int)s.size());
    bool ok=(s=="AM1 EXTERNAL GEOM");
    std::printf("\"%s\" %s\n",s.c_str(),ok?"PASS":"FAIL");
    return ok?0:1;
}
