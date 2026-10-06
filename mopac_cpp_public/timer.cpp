// timer.cpp
#include "timer.h"
#include <cstdio>
double CPU_0 = 0.0, CPU_1 = 0.0;
double wall_clock_0 = 0.0, wall_clock_1 = 0.0;
extern double second(int);
void timer(const std::string& a) {
    static bool first = true;
    static double CPU_start, CPU_last;
    static double wall_clock_last;
    if (first) {
        CPU_start = CPU_0;
        CPU_last = CPU_start;
        wall_clock_last = wall_clock_1;
        wall_clock_0 = wall_clock_last;
        first = false;
    }
    int arg1 = 1;
    double dummy = second(arg1);
    if (dummy < -300.0) return;
    double CPU_now = CPU_1;
    std::printf("%s CPU INTERVAL: %.2f, INTEGRAL: %.2f WALL-CLOCK INTERVAL: %.2f, INTEGRAL: %.2f\n",
        a.c_str(), CPU_now - CPU_last, CPU_now - CPU_start,
        wall_clock_1 - wall_clock_last, wall_clock_1 - wall_clock_0);
    CPU_last = CPU_now;
    wall_clock_last = wall_clock_1;
}
