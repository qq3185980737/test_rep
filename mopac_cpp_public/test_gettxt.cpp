// test_gettxt.cpp
#include <cstdio>
#include <string>
#include "gettxt.h"
int main() {
    std::string s = "GEO_REF=\"abc.mop\" end";
    std::string r = get_text(s, 8, 1);
    bool ok = (r == "abc.mop");
    std::printf("get_text=[%s] %s\n", r.c_str(), ok?"PASS":"FAIL");
    return ok?0:1;
}
