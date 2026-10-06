// test_win.cpp
#include <cstdio>
#include "MOPAC_for_Windows.h"
int main() { int r=mopac_win_main(0,nullptr); std::printf("win_main r=%d PASS\n",r); return 0; }
