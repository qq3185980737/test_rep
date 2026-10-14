// test_second.cpp
#include <cstdio>
#include <thread>
#include <chrono>
#include "second.h"
int main() {
    double t0=second(1);
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    double t1=second(2);
    bool ok=(t1-t0>0.04 && t1-t0<1.0);
    std::printf("dt=%g %s\n",t1-t0,ok?"PASS":"FAIL");
    return ok?0:1;
}
