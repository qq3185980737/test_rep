// test_timer.cpp
#include <cstdio>
#include "timer.h"
extern "C" double second_(int* p){ return 100.0 + *p; }
int main() { timer("step"); timer("step2"); std::printf("timer PASS\n"); return 0; }
