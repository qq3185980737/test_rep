// second.cpp
#include "second.h"
#include <chrono>
static std::chrono::steady_clock::time_point t0;
static bool initialized=false;
double second(int) {
    if (!initialized) { t0=std::chrono::steady_clock::now(); initialized=true; }
    auto t1=std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1-t0).count();
}
