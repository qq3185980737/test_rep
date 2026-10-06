// test_wt.cpp
#include <cstdio>
#include "write_trajectory.h"
int main() { double xyz[3],ch[1]; write_trajectory(xyz,1,ch,0,0,0,0); std::printf("wt links OK PASS\n"); return 0; }
