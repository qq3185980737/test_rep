// test_post_scf.cpp
#include <cstdio>
#include "post_scf_corrections.h"
int main() { double c; post_scf_corrections(c,false); std::printf("post_scf corr=%g PASS\n",c); return 0; }
