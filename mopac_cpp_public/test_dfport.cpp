// test_dfport.cpp
#include <cstdio>
#include "dfport.h"
int main() { auto s=ifport::jdate(); std::printf("jdate=%s %s\n",s.c_str(),s=="n"?"PASS":"FAIL"); return s=="n"?0:1; }
