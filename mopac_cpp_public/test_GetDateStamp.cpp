// test_GetDateStamp.cpp
#include <cstdio>
#include <string>
#include "GetDateStamp.h"
int main() { std::string d,v; GetDateStamp(d,v); bool ok=(v=="16.013"); std::printf("v=%s %s\n",v.c_str(),ok?"PASS":"FAIL"); return ok?0:1; }
