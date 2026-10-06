// test_md_C.cpp — verify md_C module defaults and array usage.
#include "md_C.h"

#include <cstdio>

using namespace md_C;

int main() {
    bool ok = true;
    if (DEF_MDSTEP != 1000) { std::printf("FAIL DEF_MDSTEP\n"); ok = false; }
    if (DEF_DTMD != 0.2) { std::printf("FAIL DEF_DTMD\n"); ok = false; }
    if (DEF_TTARGET != 16000.0) { std::printf("FAIL DEF_TTARGET\n"); ok = false; }
    if (KB != 3.166811563e-6) { std::printf("FAIL KB\n"); ok = false; }
    if (md_step_max != 1000 || dt_md != 0.2 || T_target != 16000.0) { std::printf("FAIL defaults\n"); ok = false; }
    if (l_mdzero || md_stop) { std::printf("FAIL flags false\n"); ok = false; }
    if (!md_test) { std::printf("FAIL md_test true\n"); ok = false; }
    // 2-D array usage: 3 x numat
    md_vel.assign(4, std::vector<double>(4, 0.0));
    md_vel[1][2] = 1.5;
    if (md_vel[1][2] != 1.5) { std::printf("FAIL md_vel\n"); ok = false; }
    std::printf("%s\n", ok ? "ALL PASS" : "FAILED");
    return ok ? 0 : 1;
}
