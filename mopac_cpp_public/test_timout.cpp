// test_timout.cpp
#include <cstdio>
#include <string>
#include "timout.h"
double wall_clock_0=0, wall_clock_1=3725.5;
double CPU_0=0, CPU_1=1800.25;
std::string keywrd="";
int main() { timout(6); std::printf("timout PASS\n"); return 0; }
