// test_jdate.cpp
#include <cstdio>
#include <string>
#include "jdate.h"
int main() { std::string s=jdate(); bool ok=(s.size()==8); std::printf("jdate=[%s] %s\n",s.c_str(),ok?"PASS":"FAIL"); return ok?0:1; }
